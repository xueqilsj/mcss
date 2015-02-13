from mcss import Potts

import wx
from numpy import zeros, int8, mean, std, arange, array, random
from re import split
from os import getcwd
from time import sleep, time
import threading


import matplotlib
#matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

try:
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    print "ImportError: mpl_toolkits.mplot3d"
    pass
            


def randrange(n, vmin, vmax):
    return (vmax - vmin) * random.rand(n) + vmin     
class FuncPanel(wx.Panel):
    """ A static panel for Potts
    """
    def __init__(self, parent, ID, mainFrame):
        wx.Panel.__init__(self, parent, ID)
        self.mainFrame = mainFrame;
      
        self.__do_layout()
        
        self.build_frame(parent)
        
        self.reset = False;
        self.ploting = False
        self.keeping = False
        self.job = None
        self.redraw_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_redraw_timer, self.redraw_timer)        
        self.redraw_timer.Start(500)
        self.starting_time = time()#clock() takes negative when it lasts long  
        # end wxGlade             

    def __do_layout(self):
        # begin wxGlade: MyFrame.__do_layout
        return
    def build_frame(self, mainPanel):
        self.fig = Figure((5, 5), dpi=100)
        self.axes_left = Axes3D(self.fig)
    
        self.canvas = FigCanvas(mainPanel, -1, self.fig)
        # Create the navigation toolbar, tied to the canvas
        #self.toolbar = NavigationToolbar(self.canvas)
        
        # Layout with box sizers
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        #self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        #self.vbox.AddSpacer(10)        
        self.vbox.Add(self, 0, flag=wx.ALL | wx.EXPAND)
        mainPanel.SetSizer(self.vbox)
        self.vbox.Fit(self.mainFrame)
    """Event functions"""
                    
    def on_redraw_timer(self, event):
        if self.ploting:
            return
        self.ploting = True
        self.draw_figure()
        self.mainFrame.refresh_status_message("%lg seconds" % (time() - self.starting_time))
        self.ploting = False
 
    def draw_figure(self):
        """ Redraws the figure
        """
        self.axes_left.clear()
        self.axes_left.set_axis_bgcolor('black')
        n = 100
        for c, zl, zh in [('r', -50, -25), ('b', -30, -5)]:
            xs = randrange(n, 23, 32)
            ys = randrange(n, 0, 100)
            zs = randrange(n, zl, zh)
            self.axes_left.scatter(xs, ys, zs, c=c, marker='o')
        

        self.axes_left.set_title("3D off-lattice system", size=10)
        #pylab.setp(self.axes_left.get_xticklabels(), fontsize=9)
        #pylab.setp(self.axes_left.get_yticklabels(), fontsize=9)
        self.canvas.draw()  
        

    def on_leftplot(self, event):
        curchoice = self.radio_box_leftplot.GetSelection()
        if curchoice == 0:
            self.LeftPlot = "ETraj"
        elif curchoice == 1:
            self.LeftPlot = "EHist"
       
        
    """for Mainframe calling functions"""
    def message_dlg(self, titile, msg):
            dlg = wx.MessageDialog(self, msg, titile, wx.OK | wx.ICON_INFORMATION)
            dlg.ShowModal() 
            dlg.Destroy()
            
