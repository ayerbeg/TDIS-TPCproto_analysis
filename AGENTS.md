# Analysis manual maintenance

When changing active analysis scripts, configuration, commands, or plotting behaviour:
- Preserve user edits to palettes and plot definitions.
- Update `docs/TPC_manual.tex` when the instructions or semantics change.
- Regenerate the representative run1616 ROOT plots when their appearance/configuration changes.
- Run `python docs/build_manual.py` and `python docs/build_manual.py --check`.
- Carry `docs/` and the shared `config/` with each refreshed source bundle.
The builder refreshes configuration/help/palette excerpts and examples; explanatory prose still needs maintenance.
