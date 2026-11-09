"""The editable cut DB is a data contract shared with C++/ROOT.

Run through: python analysis/tpc.py check
Fixtures use temporary DB files; they do not add cuts to real detector runs.
"""
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from detector import CutDatabase, select_event_hits, select_hit_mask, raw_histogram_spec, channel_to_pad


class RunCutTests(unittest.TestCase):
    def database(self, text):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        path = Path(directory.name) / 'cuts.ini'
        path.write_text(text)
        return CutDatabase(path)

    def test_global_run_and_separate_study_scopes(self):
        db = self.database('[global]\nreject = Time >= 228\n[run 1616]\nreject = Channel in 10,12\n')
        study = self.database('[global]\nreject = Channel >= 285\n')
        rows = dict(Channel=np.array([10, 11, 335, 11]), Time=np.array([30, 30, 30, 228]))
        np.testing.assert_array_equal(db.keep(rows, 1616), [False, True, True, False])
        np.testing.assert_array_equal(db.keep(rows, 1615), [True, True, True, False])
        np.testing.assert_array_equal(select_hit_mask([1]*4, rows['Channel'], rows['Time'],
            run=1615, cuts_db=db, voltage=True, study_cuts_db=study), [True, True, False, False])

    def test_detector_cleaning_does_not_read_the_study_database(self):
        db = self.database('[global]\nreject = Time >= 228\n')
        missing = db.path.parent / 'missing_study.ini'
        np.testing.assert_array_equal(select_hit_mask([1], [335], [30], run=1616,
            cuts_db=db, study_cuts_db=missing), [True])
        with self.assertRaises(FileNotFoundError):
            select_hit_mask([1], [335], [30], run=1616, cuts_db=db,
                            voltage=True, study_cuts_db=missing)

    def test_no_hidden_duplicate_manual_cuts(self):
        db = self.database('[global]\n')
        # Replacing the database with an empty definition must not leave
        # the former hard-coded Time/Channel cuts active behind the scenes.
        from detector import select_hit_mask
        np.testing.assert_array_equal(
            select_hit_mask([1, 1], [149, 335], [229, 500], cuts_db=db, voltage=True, study_cuts_db=db),
            [True, True])

    def test_rules_or_conditions_and_and_inclusive_endpoints(self):
        db = self.database('[run 1616]\nreject = Channel in 40..44 && Time >= 90 && Time < 95\nreject = Event in 1000..1002\n')
        rows = dict(Channel=np.array([40, 44, 44, 45, 10, 10, 10]),
                    Time=np.array([90, 94, 95, 90, 30, 30, 30]),
                    Event=np.array([1, 1, 1, 1, 999, 1000, 1002]))
        np.testing.assert_array_equal(db.keep(rows, 1616), [False, False, True, True, True, False, False])
        np.testing.assert_array_equal(db.keep(rows, 1615), np.ones(7, bool))

    def test_scalar_operators_negative_and_scientific_values(self):
        values = np.array([-5, 0, 2, 3, 4])
        for op, expected in [('==', values != 3), ('!=', values == 3),
                             ('<', values >= 3), ('<=', values > 3),
                             ('>', values <= 3), ('>=', values < 3)]:
            db = self.database(f'[global]\nreject = Charge {op} 3e0\n')
            np.testing.assert_array_equal(db.keep({'Charge': values}), expected)
        db = self.database('[global]\nreject = Charge < -1.5\nreject = Charge not in 0..4\n')
        np.testing.assert_array_equal(db.keep({'Charge': values}), [False, True, True, True, True])

    def test_aliases_and_run_number_normalization(self):
        db = self.database('[run 01616]\nreject = FADC == 153\nreject = Pad == 24\n')
        np.testing.assert_array_equal(db.keep({'Channel': np.array([152, 15, 150])}, 1616), [False, False, True])

    def test_channel_cut_propagates_to_the_mapped_pad(self):
        channels = np.array([152, 152, 150, 15, 504, 0])
        pads = channel_to_pad(channels)
        masks = [self.database(f'[run 1616]\nreject = {condition}\n').keep({'Channel': channels}, 1616)
                 for condition in ('Channel == 152', 'FADC == 153', 'Pad == 0')]
        for keep in masks:
            np.testing.assert_array_equal(keep, [False, False, True, True, True, True])
            self.assertNotIn(0, pads[keep])
            self.assertIn(1, pads[keep])
        # An electronics exclusion acts on hit rows, including disconnected
        # channels, without reinterpreting the channel number as a pad ID.
        db = self.database('[run 1616]\nreject = Channel == 0\n')
        np.testing.assert_array_equal(db.keep({'Channel': channels}, 1616), [True, True, True, True, True, False])

    def test_filter_before_first_pulse_reduction(self):
        db = self.database('[run 1616]\nreject = Channel == 42 && Time == 20\n')
        hits = select_event_hits([1, 1], [42, 42], [20, 40], [100, 200], 1, run=1616, cuts_db=db)
        self.assertEqual(hits, [(42, 200, 40)])

    def test_bad_syntax_and_missing_branch_are_not_silent(self):
        for text in ['[run1616]\nreject = Channel == 1\n', '[global]\nreject = Channnel == 1\n',
                     '[global]\nreject = Time in 95..90\n', '[global]\nreject = Charge > nan\n',
                     '[global]\nreject = Channel == 1 &&\n', '[global]\n[global]\n', '[voltage]\nreject = Channel >= 285\n']:
            with self.assertRaisesRegex(ValueError, r'cuts.ini:\d+:'):
                self.database(text)
        db = self.database('[run 1616]\nreject = Pedestal > 100\n')
        with self.assertRaisesRegex(ValueError, 'Pedestal is missing'):
            db.keep({'Channel': np.array([10])}, 1616)

    def test_raw_axis_width_range_and_log_defaults(self):
        edges, logarithmic = raw_histogram_spec(np.array([-1, 395]), 'Pad')
        self.assertEqual(len(edges) - 1, 410)
        self.assertEqual((edges[0], edges[-1]), (-.5, 409.5))
        self.assertTrue(logarithmic)
        edges, logarithmic = raw_histogram_spec(np.array([0, 575]), 'Channel')
        self.assertEqual(len(edges) - 1, 576)
        self.assertEqual((edges[0], edges[-1]), (-.5, 575.5))
        self.assertTrue(logarithmic)
        edges, _ = raw_histogram_spec(np.array([17, 229]), 'Time')
        self.assertEqual((edges[0], edges[-1]), (-.5, 229.5))
        edges, _ = raw_histogram_spec(np.array([-10, 100]), 'Charge', bins=20, limits=(0, 100))
        np.testing.assert_array_equal(edges, np.linspace(0, 100, 21))

    def test_mapped_raw_pad_counts_keep_unmapped_rows_in_underflow(self):
        # Electronics are zero-based; printed fADC153 is Channel152, Pad0.
        channels = np.array([0, 152, 150, 15, 504, 575])
        pads = channel_to_pad(channels)
        np.testing.assert_array_equal(pads, [-1, 0, 1, 24, 250, 395])
        np.testing.assert_array_equal(channels, [0, 152, 150, 15, 504, 575])
        edges, _ = raw_histogram_spec(pads, 'Pad')
        counts, _ = np.histogram(pads, edges)
        self.assertEqual(counts.sum() + (pads < edges[0]).sum(), len(channels))
        self.assertEqual(counts[0], 1)
        self.assertEqual(counts[250], 1)
        self.assertEqual((pads < edges[0]).sum(), 1)


if __name__ == '__main__':
    unittest.main()
