"""Keep config/workflow_stages.yaml in sync with the repository (entry scripts exist, Snakemake rules are defined)."""
import glob
import os
import re
import unittest

import yaml

import util

REG = os.path.join(util.PIPELINE_DIR, "config", "workflow_stages.yaml")
STATUSES = {"workflow", "tool", "project_script", "manual", "planned"}


class WorkflowStagesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(REG) as fh:
            cls.reg = yaml.safe_load(fh)
        rules = set()
        for f in glob.glob(os.path.join(util.PIPELINE_DIR, "workflow", "rules", "*.smk")):
            rules |= set(re.findall(r"^rule (\w+):", open(f).read(), re.M))
        cls.rules = rules

    def test_schema_and_unique_ids(self):
        self.assertEqual(self.reg["schema"], "g2f-workflow-stages/1")
        ids = [s["id"] for s in self.reg["stages"]]
        self.assertEqual(len(ids), len(set(ids)))

    def test_fields_and_status(self):
        for s in self.reg["stages"]:
            for k in ("id", "stage", "purpose", "status", "entry", "rules", "outputs", "gate", "human_review"):
                self.assertIn(k, s, s.get("id"))
            self.assertIn(s["status"], STATUSES, s["id"])
            if s["status"] == "workflow":
                self.assertTrue(s["rules"], f"{s['id']}: a workflow stage names its rules")
            if s["status"] == "tool":
                self.assertTrue(s["entry"], f"{s['id']}: a tool stage names its script")

    def test_entry_scripts_exist(self):
        for s in self.reg["stages"]:
            for e in s["entry"]:
                self.assertTrue(os.path.exists(os.path.join(util.PIPELINE_DIR, e)), f"{s['id']}: missing {e}")

    def test_rules_defined(self):
        for s in self.reg["stages"]:
            for r in s["rules"]:
                self.assertIn(r, self.rules, f"{s['id']}: rule {r} not in workflow/rules/*.smk")


if __name__ == "__main__":
    unittest.main()
