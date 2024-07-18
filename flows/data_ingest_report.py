# Copyright (c) 2024, Plinder Development Team
# Distributed under the terms of the Apache License 2.0
"""
The filestore instance name is: plinder-data-gen.
"""

from metaflow import FlowSpec, kubernetes, environment, step, retry

import report

MOUNT = "/plinder"

K8S = dict(
    cpu=1,
    image="ghcr.io/plinder-org/plinder:v0.1.1",
    persistent_volume_claims={
        "plinder-data-gen-pvc": MOUNT,
    },
)
ENV = dict(
    vars=dict(
        PLINDER_MOUNT=MOUNT,
        PLINDER_RELEASE="2024-04",
    )
)


class PlinderDataIngestReportFlow(FlowSpec):

    @kubernetes(**K8S)
    @environment(**ENV)
    @retry
    @step
    def start(self):
        self.dfs = report.main(upload=True)
        self.next(self.end)

    @step
    def end(self):
        for name, df in self.dfs.items():
            df.to_csv(f"reports/{name}.csv", index=False)



if __name__ == '__main__':
    PlinderDataIngestReportFlow()
