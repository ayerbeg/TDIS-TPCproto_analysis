"""2D axis definitions preserve physical pads, full extents and safe sizes."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from detector import load_2d_plots, raw_2d_axes, raw_histogram_spec


class Raw2DTests(unittest.TestCase):
    def setUp(self):
        self.arrays = dict(Pad=np.array([-1, 0, 1, 409]), Event=np.array([0, 63, 64, 65535]),
                           Hit=np.array([1, 2, 3, 4]), Charge=np.array([-10, 0, 4, 4095]),
                           Time=np.array([0, 40, 228, 229]))

    def test_requested_pairs_orientation_and_event_resolution(self):
        # Definitions are intentionally user-editable: use a fixture rather
        # than freezing the user's current plot orientation/resolution.
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'plots.csv'
            path.write_text('Name,X,Y,YBinWidth,LogZ\nEvent_vs_Pad,Pad,Event,64,true\nHit_vs_Pad,Pad,Hit,,true\nCharge_vs_Pad,Pad,Charge,,true\nTime_vs_Charge,Charge,Time,,true\n')
            plots = load_2d_plots(path)
        self.assertEqual([(p['X'], p['Y']) for p in plots],
                         [('Pad', 'Event'), ('Pad', 'Hit'), ('Pad', 'Charge'), ('Charge', 'Time')])
        x, y = raw_2d_axes(self.arrays, plots[0])
        self.assertEqual((len(x)-1, x[0], x[-1]), (410, -.5, 409.5))
        self.assertEqual((len(y)-1, y[0], y[-1]), (1024, -.5, 65535.5))
        # IDs 0..63 share a bin, while Event64 begins the next one.
        np.testing.assert_array_equal(np.searchsorted(y, [0, 63, 64], side='right'), [1, 1, 2])

    def test_one_sided_bounds_and_count_override_inherited_pad_bins(self):
        plot = dict(Name='zoom', X='Pad', Y='Time', XBinWidth='2', XMax='9.5', YMin='20')
        x, y = raw_2d_axes(self.arrays, plot)
        self.assertEqual((len(x)-1, x[0], x[-1]), (5, -.5, 9.5))
        self.assertEqual(y[0], 20)
        plot.update(XBins='2')  # Explicit count takes precedence over width.
        x, _ = raw_2d_axes(self.arrays, plot)
        self.assertEqual(len(x)-1, 2)

    def test_no_silent_massive_event_pad_allocation(self):
        with self.assertRaisesRegex(ValueError, '4 million'):
            raw_2d_axes(self.arrays, dict(Name='large', X='Pad', Y='Event', YBins='10000'))
        with self.assertRaisesRegex(ValueError, 'missing branch'):
            raw_2d_axes(self.arrays, dict(Name='bad', X='Pedestal', Y='Time'))

    def test_bad_or_duplicate_plot_definitions_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'pairs.csv'
            for rows in [[['bad/name','Pad','Time','true']],
                         [['good','Pad','Time','yes']],
                         [['good','Unknown','Time','true']],
                         [['same','Pad','Time','true'],['same','Pad','Charge','true']]]:
                with path.open('w', newline='') as stream:
                    writer = csv.writer(stream); writer.writerow(['Name','X','Y','LogZ']); writer.writerows(rows)
                with self.assertRaisesRegex(ValueError, 'pairs.csv:'):
                    load_2d_plots(path)

    def test_whitespace_around_plot_headers_values_and_empty_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'plots.csv'
            path.write_text(' Name ,\t X , Y , XBins , YBinWidth , LogZ \n test ,\t Pad , Time , 4 ,\t , true \n')
            plots = load_2d_plots(path)
            self.assertEqual((plots[0]['X'], plots[0]['Y'], plots[0]['YBinWidth']), ('Pad','Time',''))
            x, _ = raw_2d_axes(self.arrays, plots[0])
            self.assertEqual(len(x)-1, 4)
            bins = Path(directory) / 'bins.csv'
            bins.write_text(' Branch , BinWidth , Bins , Min , Max , LogY \n Pad ,\t 1 , 4 , -0.5 , 3.5 , true \n')
            edges, scale = raw_histogram_spec(np.array([0,3]), 'Pad', binning_file=bins)
            np.testing.assert_array_equal(edges, [-.5,.5,1.5,2.5,3.5])
            self.assertTrue(scale)


if __name__ == '__main__':
    unittest.main()
