"""
One-file ChEMBL target-level summary script.

It does this:

1. Read every network comparison CSV from step 02.
2. Count each ChEMBL target once within each disease.
3. Keep the contributing UniProt proteins and source row numbers.
4. Show whether any component from each target is in the disease network.
5. Save two detailed CSV files and one Excel workbook.

Input:
- results/02_network_comparisons/network_comparison_*.csv

Output:
- results/03_chembl_target_summary/Unique ChEMBL targets.csv
- results/03_chembl_target_summary/UniProt components.csv
- results/03_chembl_target_summary/chembl_target_summary.xlsx

Run:  python3 03_summarize_the_chembl_targets.py
"""

import ast
from pathlib import Path

import pandas as pd
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo


# =============================================================================
# USER INPUT
# =============================================================================

HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"
NETWORK_COMPARISON_DIR = RESULTS_DIR / "02_network_comparisons"
SUMMARY_DIR = RESULTS_DIR / "03_chembl_target_summary"

# Set this to False if you only want the four CSV files.
WRITE_EXCEL_WORKBOOK = True


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

IN_NETWORK_LEVELS = {"Tclin", "Tchem", "Tbio", "Tdark"}


def load_minimum_phase():
    """Read MINIMUM_PHASE from step 01 without running that script."""

    step_one_path = HERE / "01_hunt_the_chembl_targets.py"
    tree = ast.parse(step_one_path.read_text())

    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue

        if any(
            isinstance(target, ast.Name) and target.id == "MINIMUM_PHASE"
            for target in node.targets
        ):
            value = ast.literal_eval(node.value)
            if isinstance(value, (int, float)):
                return value

    raise SystemExit(f"Could not find MINIMUM_PHASE in {step_one_path.name}")


def phase_description(minimum_phase):
    """Return a short human-readable description of the phase filter."""

    if minimum_phase == 4:
        return "Approved drugs only"
    if minimum_phase <= 0:
        return "All recorded phases"
    return f"Phase {minimum_phase} and above"


def disease_label(disease):
    """Convert colorectal_cancer into Colorectal Cancer."""

    return disease.replace("_", " ").title()


def split_values(value):
    """Split a semicolon-separated cell into individual values."""

    if pd.isna(value) or not str(value).strip():
        return []

    return [item.strip() for item in str(value).split(";") if item.strip()]


def unique_values(values):
    """Return sorted unique non-empty values."""

    cleaned = set()

    for value in values:
        if pd.isna(value):
            continue

        value = str(value).strip()
        if value:
            cleaned.add(value)

    return sorted(cleaned)


def join_unique(values):
    """Join sorted unique values with semicolons."""

    return "; ".join(unique_values(values))


def collect_semicolon_values(series):
    """Collect unique values from a column containing semicolon-separated cells."""

    values = []

    for cell in series:
        values.extend(split_values(cell))

    return unique_values(values)


def compress_row_numbers(numbers):
    """Convert [2, 3, 4, 7] into '2-4, 7'."""

    numbers = sorted({int(number) for number in numbers})

    if not numbers:
        return ""

    ranges = []
    start = previous = numbers[0]

    for number in numbers[1:] + [None]:
        if number is not None and number == previous + 1:
            previous = number
            continue

        ranges.append(str(start) if start == previous else f"{start}-{previous}")

        if number is not None:
            start = previous = number

    return ", ".join(ranges)


def target_network_summary(group):
    """Summarize network membership for all components of one ChEMBL target."""

    levels = group["Network level"].fillna("")
    in_network = levels.isin(IN_NETWORK_LEVELS)
    not_in_network = levels.eq("not in network")
    not_mapped = levels.eq("not mapped to STRING")

    if len(group) > 0 and in_network.all():
        status = "all components in network"
    elif in_network.any():
        status = "some components in network"
    elif len(group) > 0 and not_mapped.all():
        status = "not mapped to STRING"
    else:
        status = "no components in network"

    return {
        "components_in_network": int(in_network.sum()),
        "components_not_in_network": int(not_in_network.sum()),
        "components_not_mapped": int(not_mapped.sum()),
        "target_network_status": status,
        "uniprot_components_in_network": join_unique(group.loc[in_network, "uniprot"]),
    }


def add_excel_table(sheet, start_row, end_row, end_column, table_name):
    """Add filtering and banded rows to one Excel range."""

    reference = f"A{start_row}:{get_column_letter(end_column)}{end_row}"
    table = Table(displayName=table_name, ref=reference)
    table.tableStyleInfo = TableStyleInfo(
        name="TableStyleMedium2",
        showFirstColumn=False,
        showLastColumn=False,
        showRowStripes=True,
        showColumnStripes=False,
    )
    sheet.add_table(table)


def set_column_widths(sheet, widths):
    """Apply readable widths to the columns in one worksheet."""

    for column_number, width in enumerate(widths, start=1):
        sheet.column_dimensions[get_column_letter(column_number)].width = width


def style_header_row(sheet, row_number, end_column):
    """Apply the shared dark-blue header style."""

    fill = PatternFill("solid", fgColor="1F4E78")
    border = Border(bottom=Side(style="thin", color="FFFFFF"))

    for cell in sheet[row_number][:end_column]:
        cell.fill = fill
        cell.font = Font(name="Arial", size=10, bold=True, color="FFFFFF")
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = border

    sheet.row_dimensions[row_number].height = 36


def style_data_sheet(sheet, widths, freeze_panes, table_name):
    """Apply shared formatting to a flat data worksheet."""

    sheet.sheet_view.showGridLines = False
    sheet.freeze_panes = freeze_panes
    style_header_row(sheet, 1, sheet.max_column)
    set_column_widths(sheet, widths)
    add_excel_table(sheet, 1, sheet.max_row, sheet.max_column, table_name)

    for row in sheet.iter_rows(min_row=2):
        for cell in row:
            cell.font = Font(name="Arial", size=10, color="1F2937")
            cell.alignment = Alignment(vertical="center")


# =============================================================================
# LOAD STEP 02 RESULTS
# =============================================================================


def load_component_rows():
    """Load all seven step-02 files and preserve their source row numbers."""

    paths = sorted(NETWORK_COMPARISON_DIR.glob("network_comparison_*.csv"))

    if not paths:
        raise SystemExit(
            f"No network_comparison_*.csv files found in {NETWORK_COMPARISON_DIR}. "
            "Run 02_map_targets_to_the_network.py first."
        )

    required_columns = {
        "uniprot",
        "ontology_ids",
        "chembl_target_id",
        "chembl_target",
        "chembl_target_type",
        "target_component_count",
        "action",
        "max_phase",
        "n_drugs",
        "Network level",
        "Approved drugs",
    }
    source_files = []
    frames = []

    for path in paths:
        disease = path.stem.removeprefix("network_comparison_")
        frame = pd.read_csv(path, sep=";")
        missing_columns = sorted(required_columns - set(frame.columns))

        if missing_columns:
            raise SystemExit(
                f"{path.name} is missing columns: {', '.join(missing_columns)}"
            )

        relative_path = path.relative_to(HERE).as_posix()
        source_files.append({
            "disease": disease,
            "source_file": relative_path,
            "component_rows": len(frame),
        })

        if frame.empty:
            continue

        frame.insert(0, "disease", disease)
        frame.insert(1, "source_file", relative_path)
        frame.insert(2, "source_row", frame.index + 2)
        frames.append(frame)

    components = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return components, pd.DataFrame(source_files)


# =============================================================================
# SUMMARIZE BY CHEMBL TARGET
# =============================================================================


def build_targets_by_disease(components):
    """Create one row per disease and ChEMBL target ID."""

    summary_rows = []

    for (disease, target_id), group in components.groupby(
        ["disease", "chembl_target_id"], sort=True
    ):
        drugs = collect_semicolon_values(group["Approved drugs"])
        ontology_ids = collect_semicolon_values(group["ontology_ids"])
        network = target_network_summary(group)

        summary_rows.append({
            "disease": disease_label(disease),
            "disease_key": disease,
            "chembl_target_id": target_id,
            "target_network_status": network["target_network_status"],
            "uniprot_components_in_network": network[
                "uniprot_components_in_network"
            ],
            "network_levels": join_unique(
                level
                for level in group["Network level"]
                if level in IN_NETWORK_LEVELS
            ),
            "chembl_target": join_unique(group["chembl_target"]),
            "chembl_target_type": join_unique(group["chembl_target_type"]),
            "component_count": group["uniprot"].nunique(),
            "uniprot_components": join_unique(group["uniprot"]),
            "components_in_network": network["components_in_network"],
            "components_not_in_network": network["components_not_in_network"],
            "components_not_mapped": network["components_not_mapped"],
            "action": join_unique(group["action"]),
            "max_phase": pd.to_numeric(group["max_phase"], errors="coerce").max(),
            "n_approved_drugs": len(drugs),
            "approved_drugs": "; ".join(drugs),
            "ontology_ids": "; ".join(ontology_ids),
            "source_file": group["source_file"].iloc[0],
            "source_rows": compress_row_numbers(group["source_row"]),
        })

    return pd.DataFrame(summary_rows).sort_values(
        ["disease", "chembl_target_id"], ignore_index=True
    )


def build_disease_summary(targets_by_disease, source_files):
    """Count target-level network results for every disease dataset."""

    rows = []

    for source in source_files.itertuples():
        disease_targets = targets_by_disease[
            targets_by_disease["disease_key"] == source.disease
        ]

        rows.append({
            "disease": disease_label(source.disease),
            "chembl_targets": len(disease_targets),
            "targets_with_any_component_in_network": int(
                disease_targets["components_in_network"].gt(0).sum()
            ),
            "targets_with_no_components_in_network": int(
                (
                    disease_targets["components_in_network"].eq(0)
                    & disease_targets["components_not_mapped"].lt(
                        disease_targets["component_count"]
                    )
                ).sum()
            ),
            "targets_not_mapped_to_string": int(
                (
                    disease_targets["component_count"].gt(0)
                    & disease_targets["components_not_mapped"].eq(
                        disease_targets["component_count"]
                    )
                ).sum()
            ),
            "component_rows": source.component_rows,
            "single_component_targets": int(
                disease_targets["component_count"].eq(1).sum()
            ),
            "multi_component_targets": int(
                disease_targets["component_count"].gt(1).sum()
            ),
            "source_file": source.source_file,
        })

    return pd.DataFrame(rows).sort_values("disease", ignore_index=True)


# =============================================================================
# SAVE OUTPUTS
# =============================================================================


def prepare_output_tables(targets_by_disease, components):
    """Create the two tables shared by the CSV and Excel outputs."""

    target_output = targets_by_disease[[
        "disease",
        "chembl_target_id",
        "target_network_status",
        "uniprot_components_in_network",
        "uniprot_components",
        "network_levels",
        "chembl_target",
        "chembl_target_type",
        "component_count",
        "components_in_network",
        "components_not_in_network",
        "action",
        "max_phase",
        "n_approved_drugs",
        "approved_drugs",
        "ontology_ids",
        "source_file",
        "source_rows",
    ]].rename(columns={
        "disease": "Disease",
        "chembl_target_id": "ChEMBL target ID",
        "target_network_status": "Target network status",
        "uniprot_components_in_network": "UniProt components in network",
        "network_levels": "Network levels",
        "chembl_target": "ChEMBL target",
        "chembl_target_type": "Target type",
        "component_count": "Component count",
        "uniprot_components": "UniProt components",
        "components_in_network": "Components in network",
        "components_not_in_network": "Components not in network",
        "action": "Action",
        "max_phase": "Max phase",
        "n_approved_drugs": "Approved drug count",
        "approved_drugs": "Approved drugs",
        "ontology_ids": "Ontology IDs",
        "source_file": "Source file",
        "source_rows": "Source rows",
    })
    target_status = targets_by_disease[[
        "disease_key",
        "chembl_target_id",
        "target_network_status",
    ]].rename(columns={"disease_key": "disease"})
    source_output = components.merge(
        target_status,
        on=["disease", "chembl_target_id"],
        how="left",
        validate="many_to_one",
    )
    source_output["disease"] = source_output["disease"].map(disease_label)
    source_output = source_output[[
        "disease",
        "uniprot",
        "target_network_status",
        "Network level",
        "chembl_target_id",
        "chembl_target",
        "chembl_target_type",
        "target_component_count",
        "action",
        "max_phase",
        "n_drugs",
        "Approved drugs",
        "ontology_ids",
        "source_file",
        "source_row",
    ]].rename(columns={
        "disease": "Disease",
        "uniprot": "UniProt ID",
        "target_network_status": "Target network status",
        "ontology_ids": "Ontology IDs",
        "chembl_target_id": "ChEMBL target ID",
        "chembl_target": "ChEMBL target",
        "chembl_target_type": "Target type",
        "target_component_count": "Reported component count",
        "action": "Action",
        "max_phase": "Max phase",
        "n_drugs": "Drug count",
        "source_file": "Source file",
        "source_row": "Source row",
    })

    return target_output, source_output


def write_csv_files(target_output, source_output):
    """Save the two CSV files that match the detailed Excel worksheets."""

    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    for obsolete_name in [
        "disease_summary.csv",
        "chembl_targets_by_disease.csv",
        "source_components.csv",
        "unique_chembl_targets.csv",
        "uniprot_components.csv",
    ]:
        (SUMMARY_DIR / obsolete_name).unlink(missing_ok=True)

    target_output.to_csv(SUMMARY_DIR / "Unique ChEMBL targets.csv", index=False)
    source_output.to_csv(SUMMARY_DIR / "UniProt components.csv", index=False)


def write_excel_workbook(
    disease_summary,
    target_output,
    source_output,
    minimum_phase,
):
    """Save a readable Excel workbook with the network result emphasized."""

    output_path = SUMMARY_DIR / "chembl_target_summary.xlsx"

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        workbook = writer.book

        disease_summary.to_excel(
            writer,
            sheet_name="Overview",
            index=False,
            startrow=11,
        )
        target_output.to_excel(writer, sheet_name="Unique ChEMBL targets", index=False)
        source_output.to_excel(writer, sheet_name="UniProt components", index=False)

        overview = workbook["Overview"]
        overview.sheet_view.showGridLines = False
        overview["A2"] = "ChEMBL target-level network summary"
        overview["A2"].font = Font(name="Arial", size=16, bold=True, color="1F4E78")
        for cell in overview[3][:9]:
            cell.border = Border(bottom=Side(style="medium", color="1F4E78"))

        total_targets = int(disease_summary["chembl_targets"].sum())
        targets_in = int(
            disease_summary["targets_with_any_component_in_network"].sum()
        )
        overview["A5"] = "ChEMBL targets"
        overview["B5"] = total_targets
        overview["A6"] = "Component rows"
        overview["B6"] = int(disease_summary["component_rows"].sum())
        overview["D5"] = "Targets with any component in network"
        overview["E5"] = targets_in
        overview["D6"] = "Share of targets"
        overview["E6"] = targets_in / total_targets if total_targets else 0
        overview["G5"] = "Minimum phase"
        overview["H5"] = minimum_phase
        overview["G6"] = "Included drugs"
        overview["H6"] = phase_description(minimum_phase)
        for label_cell in ["A5", "A6", "D5", "D6", "G5", "G6"]:
            overview[label_cell].font = Font(name="Arial", size=10, bold=True)

        for label_cell in ["A5", "A6"]:
            overview[label_cell].fill = PatternFill("solid", fgColor="D9EAF7")
        for label_cell in ["D5", "D6"]:
            overview[label_cell].fill = PatternFill("solid", fgColor="E2F0D9")
        for label_cell in ["G5", "G6"]:
            overview[label_cell].fill = PatternFill("solid", fgColor="FFF2CC")
        for label_cell in ["D5", "D6"]:
            overview[label_cell].alignment = Alignment(wrap_text=True, vertical="center")

        overview.row_dimensions[5].height = 30
        overview.row_dimensions[6].height = 30

        for value_cell in ["B5", "B6", "E5", "E6"]:
            overview[value_cell].font = Font(
                name="Arial", size=12, bold=True, color="1F4E78"
            )
            overview[value_cell].alignment = Alignment(horizontal="right")

        overview["H5"].font = Font(
            name="Arial", size=12, bold=True, color="1F4E78"
        )
        overview["H5"].alignment = Alignment(horizontal="right")
        overview["H6"].font = Font(
            name="Arial", size=10, bold=True, color="1F4E78"
        )

        overview["E6"].number_format = "0.0%"
        overview["A8"] = (
            "Network rule: A target is counted as in the network when at least "
            "one of its UniProt components is Tclin, Tchem, Tbio or Tdark."
        )
        overview["A9"] = (
            "Traceability: Every target retains its source CSV, original row "
            "numbers and contributing UniProt components."
        )

        for cell in [overview["A8"], overview["A9"]]:
            cell.font = Font(name="Arial", size=10, italic=True, color="1F2937")

        overview_order = [
            "disease",
            "chembl_targets",
            "targets_with_any_component_in_network",
            "targets_with_no_components_in_network",
            "multi_component_targets",
            "component_rows",
            "single_component_targets",
            "targets_not_mapped_to_string",
            "source_file",
        ]
        focused_summary = disease_summary[overview_order]

        overview_headers = [
            "Disease",
            "ChEMBL targets",
            "Targets with any component in network",
            "Targets with no components in network",
            "Multi-component targets",
            "Component rows",
            "Single-component targets",
            "Targets not mapped to STRING",
            "Source file",
        ]
        for column_number, column_name in enumerate(overview_headers, start=1):
            overview.cell(12, column_number, column_name)

        for row_number, values in enumerate(
            focused_summary.itertuples(index=False, name=None),
            start=13,
        ):
            for column_number, value in enumerate(values, start=1):
                overview.cell(row_number, column_number, value)

        style_header_row(overview, 12, len(overview_order))
        add_excel_table(
            overview,
            12,
            12 + len(focused_summary),
            len(overview_order),
            "DiseaseSummaryTable",
        )
        set_column_widths(overview, [23, 16, 24, 25, 22, 16, 22, 21, 56])

        green_fill = PatternFill("solid", fgColor="E2F0D9")
        orange_fill = PatternFill("solid", fgColor="FCE4D6")
        for row_number in range(13, 13 + len(focused_summary)):
            overview.cell(row_number, 3).fill = green_fill
            overview.cell(row_number, 4).fill = orange_fill

        overview.cell(12, 3).fill = PatternFill("solid", fgColor="548235")
        overview.cell(12, 4).fill = PatternFill("solid", fgColor="C65911")

        overview.freeze_panes = "A12"

        target_sheet = workbook["Unique ChEMBL targets"]
        style_data_sheet(
            target_sheet,
            [22, 18, 28, 38, 46, 20, 34, 24, 14, 16, 18, 18, 11, 17, 46, 38, 52, 18],
            "C2",
            "TargetsByDiseaseTable",
        )
        status_styles = {
            "all components in network": ("E2F0D9", "375623"),
            "some components in network": ("FFF2CC", "7F6000"),
            "no components in network": ("FCE4D6", "9C0006"),
            "not mapped to STRING": ("D9EAF7", "1F4E78"),
        }
        for row_number in range(2, target_sheet.max_row + 1):
            status_cell = target_sheet.cell(row_number, 3)
            colors = status_styles.get(status_cell.value)
            if colors:
                status_cell.fill = PatternFill("solid", fgColor=colors[0])
                status_cell.font = Font(name="Arial", size=10, color=colors[1])
        style_data_sheet(
            workbook["UniProt components"],
            [22, 15, 28, 20, 18, 38, 24, 20, 18, 11, 12, 52, 38, 52, 12],
            "D2",
            "SourceComponentsTable",
        )
        component_sheet = workbook["UniProt components"]
        for row_number in range(2, component_sheet.max_row + 1):
            status_cell = component_sheet.cell(row_number, 3)
            colors = status_styles.get(status_cell.value)
            if colors:
                status_cell.fill = PatternFill("solid", fgColor=colors[0])
                status_cell.font = Font(name="Arial", size=10, color=colors[1])

    return output_path


# =============================================================================
# RUN EVERYTHING
# =============================================================================


def main():
    minimum_phase = load_minimum_phase()
    components, source_files = load_component_rows()
    targets_by_disease = build_targets_by_disease(components)
    disease_summary = build_disease_summary(targets_by_disease, source_files)
    target_output, source_output = prepare_output_tables(
        targets_by_disease,
        components,
    )

    write_csv_files(target_output, source_output)

    workbook_path = None
    if WRITE_EXCEL_WORKBOOK:
        workbook_path = write_excel_workbook(
            disease_summary,
            target_output,
            source_output,
            minimum_phase,
        )

    print(
        f"Minimum phase:              {minimum_phase} "
        f"({phase_description(minimum_phase)})"
    )
    print(f"Component rows:              {len(components)}")
    print(f"Disease-target records:      {len(targets_by_disease)}")
    print("Saved CSV files:             Unique ChEMBL targets.csv, "
          "UniProt components.csv")
    print(f"Output folder:               {SUMMARY_DIR}")

    if workbook_path:
        print(f"Saved Excel workbook:        {workbook_path.name}")


if __name__ == "__main__":
    main()
