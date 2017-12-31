from click.testing import CliRunner

import os
import pkg_resources
import tempfile
import shutil
import unittest

import pyclone.cli


class Test(unittest.TestCase):

    def setUp(self):
        self.data_file = pkg_resources.resource_filename('pyclone', 'test/data/tiny.tsv')

        self.tmp_dir = tempfile.mkdtemp()

        self.trace_file = os.path.join(self.tmp_dir, 'trace.h5')

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_instantiated_new(self):
        self._run_instantiated_sampler()

        self._test_file_exists(self.trace_file)

    def test_instantiated_plot(self):
        self._run_instantiated_sampler()

        self._test_plot_commands()

    def test_instantiated_table(self):
        self._run_instantiated_sampler()

        self._test_table_commands()

    def test_marginal_new(self):
        self._run_marginal_sampler()

        self._test_file_exists(self.trace_file)

    def test_marginal_plot(self):
        self._run_marginal_sampler()

        self._test_plot_commands()

    def test_marginal_table(self):
        self._run_marginal_sampler()

        self._test_table_commands()

    def _run_instantiated_sampler(self):
        cmd = ['analysis', 'new', '-i', self.data_file, '-t', self.trace_file, '-n', 1]

        self._run_pyclone_cmd(cmd)

    def _run_marginal_sampler(self):
        cmd = ['analysis', 'new', '-i', self.data_file, '-t', self.trace_file, '-n', 1, '--grid-size', 10]

        self._run_pyclone_cmd(cmd)

    def _run_pyclone_cmd(self, cmd):
        runner = CliRunner()

        result = runner.invoke(pyclone.cli.pyclone, cmd)

        try:
            self.assertEqual(result.exit_code, 0)

        except AssertionError as e:
            print(cmd)

            print(result.output)

            raise e

    def _test_file_exists(self, path):
        self.assertTrue(os.path.exists(path))

    def _test_plot_commands(self):
        out_file = os.path.join(self.tmp_dir, 'test.png')

        for plot_type in ['density', 'line', 'scatter']:
            cmd = [
                'post-process', 'plot', 'clusters',
                '-t', self.trace_file, '-o', out_file, '-f', plot_type, '--grid-size', 10
            ]

            self._run_pyclone_cmd(cmd)

            self._test_file_exists(out_file)

            os.unlink(out_file)

        for plot_type in ['ccf-density', 'ccf-line', 'ccf-scatter', 'similarity-matrix', 'vaf-line', 'vaf-scatter']:
            cmd = [
                'post-process', 'plot', 'loci',
                '-t', self.trace_file, '-o', out_file, '-f', plot_type
            ]

            self._run_pyclone_cmd(cmd)

            self._test_file_exists(out_file)

            os.unlink(out_file)

    def _test_table_commands(self):
        out_file = os.path.join(self.tmp_dir, 'test.tsv')

        for table_type in ['cluster', 'loci', 'old']:
            cmd = [
                'post-process', 'table',
                '-t', self.trace_file, '-o', out_file, '-f', table_type, '--grid-size', 10
            ]

            self._run_pyclone_cmd(cmd)

            self._test_file_exists(out_file)

            os.unlink(out_file)


if __name__ == "__main__":
    unittest.main()
