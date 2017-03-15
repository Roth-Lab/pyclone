from pyclone.sampler import DataPoint, DirichletProcessSampler, PartitionSampler
from pyclone.trace import TraceDB

import numpy as np
import shutil


def main():
    num_data_points = 8

    num_iters = 1000

    out_dir = 'tmp'

    tumour_content = 1.0

    alpha_shape = 1

    alpha_rate = 1

    eps = 1e-3

    data_points = {}

    for i in range(num_data_points):
        cn_r = np.array([2, 2])
        cn_v = np.array([2, 2])

        mu_v = np.array([0.5, 0.999])

        log_pi = np.array([np.log(0.5), np.log(0.5)])

        if i <= 4:
            b = 1000
            d = 2000
        else:
            b = 2000
            d = 2000

        data_points[str(i)] = DataPoint(b, d, eps, cn_r, cn_v, mu_v, log_pi)

    results_db = TraceDB(out_dir, data_points.keys())

    sampler = DirichletProcessSampler(tumour_content, None, alpha_shape, alpha_rate)

    sampler.sample(data_points.values(), results_db, num_iters, 100, seed=1)

    shutil.rmtree(out_dir)

if __name__ == "__main__":
    import line_profiler

    profiler = line_profiler.LineProfiler(DirichletProcessSampler.interactive_sample,
                                          PartitionSampler.sample,
                                          DataPoint.log_p)

    profiler.run("main()")

    profiler.print_stats()

#    import cProfile
#
#    cProfile.run("main()")

#    main()
