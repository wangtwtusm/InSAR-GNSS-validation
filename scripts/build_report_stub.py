#!/usr/bin/env python3
"""Build a metrics-grounded Markdown starting point for a validation report.

The generated document deliberately avoids inventing scientific conclusions.
Use references/report-guide.md to turn the evidence bundle into a reviewed
technical report or reviewer response.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


def parse_labeled_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Use LABEL=PATH for --metric-table")
    label, path_text = value.split("=", 1)
    if not label.strip() or not path_text.strip():
        raise argparse.ArgumentTypeError("Both LABEL and PATH are required")
    return label.strip(), Path(path_text).expanduser()


def markdown_escape(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def csv_to_markdown(path: Path) -> str:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return "_No rows._"
    columns = list(rows[0])
    lines = [
        "| " + " | ".join(markdown_escape(column) for column in columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append(
            "| " + " | ".join(markdown_escape(row.get(column, "")) for column in columns) + " |"
        )
    return "\n".join(lines)


def metadata_table(metadata: dict[str, Any]) -> str:
    preferred = [
        "comparison_role",
        "reference_frame",
        "observation_window",
        "velocity_units",
        "vertical_positive",
        "coordinate_crs",
        "insar_product_stage",
        "los_sign",
        "incidence_definition",
        "heading_convention",
        "angle_units",
    ]
    rows = [(key, metadata[key]) for key in preferred if key in metadata]
    rows.extend(
        (key, metadata[key])
        for key in sorted(metadata)
        if key not in {name for name, _ in rows}
    )
    return "\n".join(
        ["| Field | Declared value |", "| --- | --- |"]
        + [f"| {markdown_escape(key)} | {markdown_escape(value)} |" for key, value in rows]
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a neutral Markdown evidence bundle from saved validation metrics."
    )
    parser.add_argument("--title", required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument(
        "--metric-table",
        action="append",
        type=parse_labeled_path,
        default=[],
        metavar="LABEL=PATH",
        help="Repeat for each saved CSV table to include.",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not args.metadata.is_file():
        raise FileNotFoundError(args.metadata)
    for _, path in args.metric_table:
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite: {args.output}")

    with args.metadata.open(encoding="utf-8") as stream:
        metadata = json.load(stream)

    lines = [
        f"# {args.title}",
        "",
        "## Methodological and product definition",
        "",
        metadata_table(metadata),
        "",
        "## Executive evidence summary",
        "",
        "The tables below preserve datum-retaining bias and raw RMS. Centered RMS, when present, is a secondary spatial-scatter diagnostic and does not replace raw accuracy. Interpret each result according to whether the comparison is internal, out of sample, or independent.",
        "",
        "## Scope, products, and statistical conventions",
        "",
        "Residuals are defined as InSAR minus GNSS. Strict statistics use station-containing finite pixels and the declared quality-control rules; unavailable or excluded stations do not enter the metric sample. No constant offset is removed unless a result is explicitly labeled as a centered diagnostic.",
        "",
        "## Validation evidence",
        "",
    ]
    if args.metric_table:
        for label, path in args.metric_table:
            lines.extend([f"### {label}", "", csv_to_markdown(path), ""])
    else:
        lines.extend(
            [
                "No metric table was supplied. Add saved numerical outputs before drawing conclusions.",
                "",
            ]
        )
    lines.extend(
        [
            "## Interpretation boundaries",
            "",
            "- State the evidence tier and sample definition with every headline result.",
            "- Do not call fitting-station agreement independent accuracy.",
            "- Do not use correlation or centered RMS to dismiss a nonzero bias.",
            "- Do not extend validation claims beyond the sampled spatial and temporal domain.",
            "- Treat reference-frame, observation-window, estimator, and equipment-history differences as possible contributors when products differ.",
            "",
            "## Data provenance and reproducibility",
            "",
            "The final report should list frozen input snapshots, hashes, product stages, selection thresholds, exclusions, software versions, and random seeds. Numerical tables should be read from saved validation outputs rather than transcribed from figures.",
            "",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
