"""Script to test and explore the issue with bad repainting of a the SplitterWindow."""
import wx

import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas

class CanvasTest(Canvas):
    def _onPaint(self, evt):
        drawDC = wx.PaintDC(self)
        if not self._isDrawn:
            self.draw(drawDC=drawDC)
        else:
            self.gui_repaint(drawDC=drawDC)
         #evt.Skip()


class Plot(wx.Panel):
    def __init__(self, parent, id=-1, dpi=None, **kwargs):
        wx.Panel.__init__(self, parent, id=id, **kwargs)
        self.figure = mpl.figure.Figure(dpi=dpi, figsize=(2, 2))
        self.canvas = Canvas(self, -1, self.figure)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        self.SetSizer(sizer)


class MainFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        wx.Frame.__init__(self, *args, **kwds)

        self.ver_splitter = wx.SplitterWindow(self, -1, style=wx.SP_3D)
        self.hor_splitter = wx.SplitterWindow(self.ver_splitter, -1, style=wx.SP_3D)

        self.input_notebook = wx.Notebook(self.hor_splitter, -1, style=wx.NB_BOTTOM | wx.NO_BORDER)
        self.plot_notebook = wx.Notebook(self.hor_splitter, -1, style=wx.NB_BOTTOM | wx.NO_BORDER)

        self.panel_left = wx.Notebook(self.ver_splitter)

        self.plot_data = Plot(self.plot_notebook)
        self.plot_notebook.AddPage(self.plot_data, "Data")


        self.ver_splitter.SplitVertically(self.panel_left, self.hor_splitter)
        self.hor_splitter.SplitHorizontally(self.plot_notebook, self.input_notebook)

        self.Layout()

        axes = self.plot_data.figure.gca()
        axes.plot([1, 2, 3], [1, 10, 100.])
        #axes.set_yticklabels(['$100$', '$10$', '$20$'])
        axes.set_yticklabels(['100', '10', '20'])
        #axes.set_yscale('log')

if __name__ == "__main__":
    app = wx.App(False)
    main_frame = MainFrame(None)
    main_frame.Show()
    app.MainLoop()