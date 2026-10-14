# config/local/ — machine and project settings (not tracked in Git)

Everything here except this README and `local.example.yaml` is untracked.

* `config.yaml` (optional): the project's top-level configuration. If present
  it is used instead of `config/config.yaml` (unless `$PIPELINE_CONFIG` names
  another file). Start from a copy of `config/config.yaml`, then set e.g.
  `tools:` (such as `config/tools/biohpc.yaml`), `resources:`, `legacy_dirs:`,
  `raw_manifest:`, `remote_copy_table:` and the `launcher:` defaults.
* `local.yaml`: copy `local.example.yaml` and edit it. It holds:
  * `external_inputs`: read-only references to read files that live outside
    `data/` (e.g. on another filesystem). They are inspected exactly like
    files in `data/<SAMPLE>/{hifi,ont}/` and are never moved or modified.
  * `busco_lineage_dir`: directory containing BUSCO lineage datasets.
  * `containers_dir`: optional override of the tools file's container path.
* Machine-specific resources files (e.g. `resources/my_server.yaml`), referred
  to from `config.yaml`.
