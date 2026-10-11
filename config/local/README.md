# config/local/ — machine-specific settings (not tracked in Git)

Copy `local.example.yaml` to `local.yaml` and edit it. `local.yaml` holds:

* `external_inputs`: read-only references to read files that live outside
  `data/` (e.g. legacy PHBW8 reads). They are inspected exactly like files in
  `data/<SAMPLE>/{hifi,ont}/` and are never moved or modified.
* `busco_lineage_dir`: directory containing BUSCO lineage datasets.
* `containers_dir`: optional override of the tools file's container path.

Only this README and `local.example.yaml` are tracked.
