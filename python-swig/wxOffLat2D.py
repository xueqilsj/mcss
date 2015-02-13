from mcss import LqFluid, PyJackknife_MeanStdErrTau_D

import wx
from numpy import zeros, int8, mean, std, arange, array, random
from re import split
from os import getcwd
from time import strftime, time
import threading


import matplotlib
#matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
import pylab

class CFluid2D(LqFluid):
    def __init__(self):
        LqFluid.__init__(self)#must use 'self'
        self.RefreshSamp = 100
        self.Stride = 1
        self.eps = 119; self.sig = 0.341;self.rs = 8.5; self.rc = 10
        self.Corxs = None; self.Corys = None; self.Fxs = None; self.Fys = None
    def Reset(self, dim0, dim1, ds, dmax):
        self.InitParameters(ds, dmax)
        self.Corxs, self.Corys, self.Fxs, self.Fys = self.PyGetCor(dim0, dim1)
        self.InitConf(1)
        self.RDF = zeros((self.nN * (self.nN - 1) / 2, 1), float) #for self.potts.PyTraj2Hist(ydata, self.potts.H)
        self.TR = []#for trajectory
        self.TRX = []#for trajectory

    def RunMetroplisTrial(self):

        for CurSample in range(self.RefreshSamp):
            self.MetroplisSweepGCE(self.Stride)
            self.TR.append(self.GetArg(self.AI_HAMI))



def randrange(n, vmin, vmax):
    return (vmax - vmin) * random.rand(n) + vmin     
class FuncPanel(wx.Panel):
    """ A static panel for Potts
    """
    def __init__(self, parent, ID, mainFrame):
        self.lqfluid = CFluid2D()
        wx.Panel.__init__(self, parent, ID)
        self.mainFrame = mainFrame;
        self.radio_regular = wx.RadioButton(self, -1, label="Regular", style=wx.RB_GROUP)
        self.label_l0l1 = wx.StaticText(self, -1, "L0,L1,dL,dMax:", style=wx.ALIGN_CENTRE)
        self.text_ctrl_l0l1 = wx.TextCtrl(self, -1, "10,10,0.47,0.05")
        self.radio_backup = wx.RadioButton(self, -1, label="Backup")
            
        
        self.label_stride = wx.StaticText(self, -1, "Stride:")
        self.text_ctrl_stride = wx.TextCtrl(self, -1, "1")
        self.label_refreshSamp = wx.StaticText(self, -1, "Ref. Samples:")
        self.text_ctrl_refreshSamp = wx.TextCtrl(self, -1, "80")
        self.label_bau = wx.StaticText(self, -1, "Beta,alpha,U/N:")
        self.text_ctrl_bau = wx.TextCtrl(self, -1, "0.25,0,-3")
        self.label_lqPar = wx.StaticText(self, -1, "eps,sig,rs,rc:")
        self.text_ctrl_lqPar = wx.TextCtrl(self, -1, "119,0.341,7.5,9")

        self.label_range = wx.StaticText(self, -1, "From:To ", style=wx.ALIGN_CENTRE)
        self.text_ctrl_range = wx.TextCtrl(self, -1, "1:")
    
        self.radio_box_upperplot = wx.RadioBox(self, -1, "Upper subplot", choices=["Energy trajectory", "Energy histogram", "Radial DF"], majorDimension=3, style=wx.RA_SPECIFY_ROWS)
        self.button_start = wx.Button(self, -1, "Start")
        self.button_analysis = wx.Button(self, -1, "Analysis")
        self.button_reset = wx.Button(self, -1, "Reset")
        
        self.__do_layout()
        
        self.build_frame(parent)
        
        self.reset = False;
        self.ploting = False
        self.keeping = False
        self.job = None
        self.record_figure_checked = False   
        self.transX = matplotlib.transforms.blended_transform_factory(self.axes_upper.transData, self.axes_upper.transAxes)    
        self.transY = matplotlib.transforms.blended_transform_factory(self.axes_upper.transAxes, self.axes_upper.transData)

        self.redraw_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_redraw_timer, self.redraw_timer)        
        self.redraw_timer.Start(300)
        self.starting_time = time()#clock() takes negative when it lasts long  
        self.radio_regular.SetFocus()
        # end wxGlade             

    def __do_layout(self):
   # begin wxGlade: MyFrame.__do_layout
              
        self.grp1_box = wx.StaticBox(self, -1, "Lennard-Jones Fluid")
        self.grp1_box.SetMinSize((220, 100))
        grp1_sizer = wx.StaticBoxSizer(self.grp1_box, wx.VERTICAL)

        grp1_sizer.Add(self.radio_regular, 0, wx.ALL, 1)
        sizer_l0l1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_l0l1.Add(self.label_l0l1, 0, wx.ALL, 3)
        sizer_l0l1.Add(self.text_ctrl_l0l1, 0, wx.ALL, 1)
        self.text_ctrl_l0l1.SetMinSize((105, 26))
        grp1_sizer.Add(sizer_l0l1, 0, wx.ALL, 1)
        grp1_sizer.Add(self.radio_backup, 0, wx.ALL, 1)
        
        self.radio_regular.Bind(wx.EVT_RADIOBUTTON, self.on_regular)
        self.text_ctrl_l0l1.Bind(wx.EVT_TEXT, self.on_l0l1)
        self.radio_backup.Bind(wx.EVT_RADIOBUTTON, self.on_backup)
        
        grp3_box = wx.StaticBox(self, -1, "")
        grp3_box.SetMinSize((200, 100))
        grp3_sizer = wx.StaticBoxSizer(grp3_box, wx.HORIZONTAL)
        
        sizer_sample1 = wx.BoxSizer(wx.VERTICAL)
        sizer_sample1.Add(self.label_stride, 0, wx.ALL, 6)
        sizer_sample1.Add(self.label_refreshSamp, 0, wx.ALL, 6)
        sizer_sample1.Add(self.label_bau, 0, wx.ALL, 6)
        sizer_sample1.Add(self.label_lqPar, 0, wx.ALL, 6)
        sizer_sample2 = wx.BoxSizer(wx.VERTICAL)
        sizer_sample2.Add(self.text_ctrl_stride, 0, wx.ALL, 2)
        sizer_sample2.Add(self.text_ctrl_refreshSamp, 0, wx.ALL, 2)
        sizer_sample2.Add(self.text_ctrl_bau, 0, wx.ALL, 2)
        sizer_sample2.Add(self.text_ctrl_lqPar, 0, wx.ALL, 2)
        self.text_ctrl_stride.SetMinSize((96, 26))
        self.text_ctrl_refreshSamp.SetMinSize((96, 26))
        self.text_ctrl_bau.SetMinSize((96, 26))
        self.text_ctrl_lqPar.SetMinSize((96, 26))
        self.text_ctrl_stride.Bind(wx.EVT_TEXT, self.on_stride)
        self.text_ctrl_refreshSamp.Bind(wx.EVT_TEXT, self.on_refreshSamp)
        self.text_ctrl_bau.Bind(wx.EVT_TEXT, self.on_bau)
        self.text_ctrl_lqPar.Bind(wx.EVT_TEXT, self.on_lqPar)
        
        grp3_sizer.Add(sizer_sample1, 0, wx.ALL, 2)
        grp3_sizer.Add(sizer_sample2, 0, wx.ALL, 2)
       
        
        grp4_sizer = wx.BoxSizer(wx.VERTICAL)
        grp4_sizer.Add(self.radio_box_upperplot, 0, wx.ALL, 2)
        self.radio_box_upperplot.SetMinSize((220, 100))
        self.radio_box_upperplot.Bind(wx.EVT_RADIOBOX, self.on_upperplot)
        self.radio_box_upperplot.SetSelection(0)
        #self.on_upperplot(None)        
        grp5_box = wx.StaticBox(self, -1, "Trajectory Range")
        grp5_sizer = wx.StaticBoxSizer(grp5_box, wx.HORIZONTAL)  
        grp5_box.SetMinSize((220, 30))
        grp5_sizer.Add(self.label_range, 0, wx.ALL, 6)
        grp5_sizer.Add(self.text_ctrl_range, 0, wx.ALL, 1)
        self.text_ctrl_range.Bind(wx.EVT_TEXT, self.on_range)
        grp4_sizer.Add(grp5_sizer, 0, wx.ALL, 2)

 
        grp6_box = wx.StaticBox(self, -1, "")
        grp6_sizer = wx.StaticBoxSizer(grp6_box, wx.HORIZONTAL)
        grp6_box.SetMinSize((220, 30))
        self.button_reset.SetMinSize((61, 30))
        self.button_analysis.SetMinSize((65, 30))
        self.button_start.SetMinSize((61, 30))
        grp6_sizer.Add(self.button_reset, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        grp6_sizer.Add(self.button_analysis, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        grp6_sizer.Add(self.button_start, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
           
        self.button_start.Bind(wx.EVT_BUTTON, self.on_start)
        self.button_reset.Bind(wx.EVT_BUTTON, self.on_reset)
        self.button_analysis.Bind(wx.EVT_BUTTON, self.on_analysis)
        self.button_start.Disable()  
        self.button_analysis.Disable()
        self.sizer_body = wx.BoxSizer(wx.VERTICAL)
        self.sizer_body.Add(grp1_sizer, border=5, flag=wx.ALL)
        self.sizer_body.Add(grp3_sizer, border=5, flag=wx.ALL)
        #self.sizer_body.AddSpacer(24)
        self.sizer_body.Add(grp4_sizer, border=5, flag=wx.ALL)
        self.sizer_body.Add(grp6_sizer, border=5, flag=wx.ALL)
        self.SetSizer(self.sizer_body)
        self.sizer_body.Fit(self)
        return
    def build_frame(self, mainPanel):
       
        self.fig = Figure((4.8, 6.5), dpi=100)
        self.axes_upper = self.fig.add_subplot(211)
        self.axes_lower = self.fig.add_subplot(212)
        self.axes_upper.set_axis_bgcolor('burlywood')
        self.axes_lower.set_axis_bgcolor('burlywood')
        self.fig.subplots_adjust(bottom=0.06, top=0.95, left=0.18, right=0.95)
        #self.fig.subplots_adjust(top=0.85)
        self.canvas = FigCanvas(mainPanel, -1, self.fig)
        # Bind the 'pick' event for clicking on one of the bars
        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)
        
        # Layout with box sizers
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
                
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        #self.vbox.AddSpacer(10)    
        self.hbox.Add(self.vbox, 0, flag=wx.ALL | wx.EXPAND)
        self.hbox.Add(self, 0, flag=wx.ALL | wx.EXPAND)
        mainPanel.SetSizer(self.hbox)
        self.hbox.Fit(self.mainFrame)
        
    def on_reset(self, event):
        self.reset = False;
        self.button_start.Disable()
        self.on_l0l1(None)
        try:
            l0 = int(self.l0l1[0])
            l1 = int(self.l0l1[1])
            dl = float(self.l0l1[2])
            dmax = float(self.l0l1[3])
        except:
            self.mainFrame.refresh_status_message("Unknown L0,L1,dL,dMax:" + ",".join(self.l0l1))
            return
        
        if l0 < 0 or l1 < 0 or dl < 0 or dmax < 0:
            self.mainFrame.refresh_status_message("Illegal L0,L1,dL,dMax:" + ",".join(self.l0l1))
            return

        self.lqfluid.Reset(l0, l1, dl, dmax)
        self.on_lqPar(None)
        self.lqfluid.Gauge()
        self.on_refreshSamp(None)
        self.on_stride(None)
        self.on_bau(None)
        self.on_range(None)
        self.on_upperplot(None)

        self.button_start.Enable()

        self.axes_upper.clear()
        self.axes_upper.set_xlabel('r', fontsize=9)
        self.axes_upper.set_ylabel(r'$U_{lq}(r)$', fontsize=9)
        self.axes_upper.set_axis_bgcolor('black')
        self.axes_upper.grid(True, color='gray')
        pylab.setp(self.axes_upper.get_xticklabels(), fontsize=9)
        pylab.setp(self.axes_upper.get_yticklabels(), fontsize=9)
        XY = self.lqfluid.PyGetXY_LQ(self.lqfluid.sigma / 3, 5 * self.lqfluid.sigma  , 80)
        self.axes_upper.plot(XY[0], XY[1], 'g-');
        self.axes_upper.text(0.05, 0.9, r"$\epsilon=%g, \sigma=%g,r_s=%g,r_c=%g$" % (self.lqfluid.epsilon, self.lqfluid.sigma, self.lqfluid.rs, self.lqfluid.rc), \
                             horizontalalignment='left', verticalalignment='top', color="red", size=12, transform=self.axes_upper.transAxes)

        self.axes_lower.clear()
        self.axes_lower.set_axis_bgcolor('black')
        self.axes_lower.set_title("Lennard-Jones fluid,N=%d" % self.lqfluid.nN, size=10)
        # for configuration plotting
        self.axes_lower.plot(self.lqfluid.Corxs, self.lqfluid.Corys, 'go')            
        #self.axes_lower.set_xlim((-self.lqfluid.pbcx / 2.0, self.lqfluid.pbcx / 2.0))
        self.canvas.draw()  
        
        self.reset = True;
        return True
                 
    def on_start(self, event):
        if not self.keeping:
            self.keeping = True;
            label = "Pause"
            self.button_analysis.Disable()
            self.button_reset.Disable()
            self.starting_time = time()#clock() takes negative when it lasts long  
        else:
            self.keeping = False;
            label = "Start"
            self.button_analysis.Enable()
            self.button_reset.Enable()
        self.button_start.SetLabel(label)
        
    
    def on_redraw_timer(self, event):
        if self.ploting:
            return
        if not self.keeping:
            return
        
        if  self.job == None:
            self.job = threading.Thread(target=self.lqfluid.RunMetroplisTrial)
            self.job.start()
        else:
            if not self.job.isAlive():
                self.ploting = True
                self.draw_figure()
                self.mainFrame.refresh_status_message("%lg seconds" % (time() - self.starting_time))
                self.job = threading.Thread(target=self.lqfluid.RunMetroplisTrial)
                self.job.start()
                self.ploting = False
    def get_data(self, TR):
        if len(self.range) == 1:
            return False
        try:
            f = int(self.range[0])
        except:
            self.mainFrame.refresh_status_message("Unknown range:" + ":".join(self.range))
            return False
        
        if self.range[1] == '':
            self.ydata = TR[f:]
            self.tdata = range(len(TR) - len(self.ydata), len(TR));
        else:
            try:
                to = int(self.range[1])
            except:
                self.mainFrame.refresh_status_message("Unknown range:" + ":".join(self.range))
                return False

            self.ydata = TR[f:to]
            x0 = len(TR) - len(self.ydata) if f < 0 else f
            self.tdata = range(x0, x0 + len(self.ydata))
        if len(self.ydata) == 0:
            return False
        return True   
    def on_analysis(self, event):
        #self.lqfluid.Gauge()
        #print self.lqfluid.TR[-1], self.lqfluid.Evdw
        if self.upperPlot == 'ETraj':
            if  len(self.ydata) == 0:
                return
            self.draw_ananlysis(r'$\bar E=%g\pm%g$, $\tau_i=%g$', self.ydata, Y=True)
        elif self.upperPlot == 'EHist':
            if  len(self.ydata) == 0:
                return
            self.draw_ananlysis(r'$\bar E=%g\pm%g$, $\tau_i=%g$', self.ydata, Y=False)
        elif self.upperPlot == 'RDF':
            """
            rdf, ni = self.lqfluid.PyGetRDF(self.lqfluid.RDF)#save to self.lqfluid.RDF
            if ni < 0:
                return
            print ni
            self.draw_ananlysis(r'$\bar r=%g\pm%g$, $\tau_i=%g$', rdf, Y=False, dlen=20)
            """
        return
    def draw_ananlysis(self, data, labels, Y, tlen=80):
            m, s, merr , ptau, taui, tauCorr = PyJackknife_MeanStdErrTau_D(data, 10, tlen)
            if taui == -1:
                return 
            print "statistic of m, s, merr , ptau, taui, tauCorr:", m, s, merr , ptau, taui
            print tauCorr
            if Y:
                rect = matplotlib.patches.Rectangle((0, m - merr), width=1, height=2 * merr, transform=self.transY, facecolor='yellow', alpha=0.5)
                l1 = matplotlib.lines.Line2D([0, 1], [m, m], transform=self.transY, linestyle="-.", color="red")
            else:
                rect = matplotlib.patches.Rectangle((m - merr, 0), width=2 * merr, height=1, transform=self.transX, facecolor='yellow', alpha=0.5)
                l1 = matplotlib.lines.Line2D([m, m], [0, 1], transform=self.transX, linestyle="-.", color="red")
  
            self.axes_upper.text(0.05, 0.9, labels % (m, merr, ptau), horizontalalignment='left', verticalalignment='top', color="red", size=10, transform=self.axes_upper.transAxes)
            self.axes_upper.lines.extend([l1])
            self.axes_upper.add_patch(rect)
            self.canvas.draw()   
            self.record_figure()
    def draw_figure1(self):
        """ Redraws the figure
        """
        self.axes_left.clear()
        n = 100
        for c, zl, zh in [('r', -50, -25), ('b', -30, -5)]:
            xs = randrange(n, 23, 32)
            ys = randrange(n, 0, 100)
            zs = randrange(n, zl, zh)
            self.axes_left.scatter(xs, ys, zs, c=c, marker='o')
            
        #pylab.setp(self.axes_left.get_xticklabels(), fontsize=9)
        #pylab.setp(self.axes_left.get_yticklabels(), fontsize=9)
        self.canvas.draw()  
        
    def draw_figure(self):
        """ Redraws the figure
        """
        
        
        self.axes_upper.clear()
        #self.axes_left.lines.remove(self.axes_left_line)
        self.axes_upper.set_title(r"LJ-fluid (GCE),$N=%d$,$\beta=%lg$,$\alpha=%lg$,$U/N=%lg$,$dMax=%.4lf$" % (self.lqfluid.nN, self.lqfluid.dBeta, self.lqfluid.dAlpha, self.lqfluid.dUn, self.lqfluid.dMax), size=10)
        self.axes_upper.grid(True, color='gray')
        self.axes_upper.set_axis_bgcolor('black')
        if self.upperPlot == 'ETraj':
            if not self.get_data(self.lqfluid.TR):
                return
            self.axes_upper.set_xlabel('Samples', fontsize=9)
            self.axes_upper.set_ylabel('Energy', fontsize=9)
            self.axes_upper.plot(self.tdata, self.ydata)
        elif self.upperPlot == 'EHist':
            if not self.get_data(self.lqfluid.TR):
                return
            self.axes_upper.set_xlabel('Energy', fontsize=9)
            self.axes_upper.set_ylabel('Histogram', fontsize=9)
            n, bins, patches = self.axes_upper.hist(self.ydata, 1000, histtype='step', normed=False)
            pylab.setp(patches, facecolor='r', alpha=0.75)
        elif self.upperPlot == 'RDF':
            rdf, ni = self.lqfluid.PyGetRDF(self.lqfluid.RDF)#save to self.lqfluid.RDF
            if ni < 0:
                return
            self.axes_upper.set_xlabel('r', fontsize=9)
            self.axes_upper.set_ylabel('Histogram', fontsize=9)
            n, bins, patches = self.axes_upper.hist(rdf, self.lqfluid.nN * 2, histtype='step', normed=False)
            pylab.setp(patches, facecolor='r', alpha=0.75)
        for label in self.axes_upper.xaxis.get_ticklabels() + self.axes_upper.get_yticklabels():
            # label is a Text instance
            #label.set_color('red')
            label.set_rotation(20)
            label.set_fontsize(9)

        #pylab.setp(self.axes_left.get_xticklabels(), fontsize=9)
        #pylab.setp(, fontsize=9)
        
        self.axes_lower.clear()
        self.axes_lower.plot(self.lqfluid.Corxs, self.lqfluid.Corys, 'go')
        self.axes_lower.set_xlabel('x', fontsize=9)
        self.axes_lower.set_ylabel('y', fontsize=9)  
        self.axes_lower.set_axis_bgcolor('black')
        self.axes_lower.set_xlim((-self.lqfluid.pbcx / 2.0, self.lqfluid.pbcx / 2.0))
        self.canvas.draw()  
        self.record_figure()  
    """for Mainframe calling functions"""
    def message_dlg(self, titile, msg):
            dlg = wx.MessageDialog(self, msg, titile, wx.OK | wx.ICON_INFORMATION)
            dlg.ShowModal() 
            dlg.Destroy()
            
    def on_upperplot(self, event):
        curchoice = self.radio_box_upperplot.GetSelection()
        if curchoice == 0:
            self.upperPlot = "ETraj"
        elif curchoice == 1:
            self.upperPlot = "EHist"
        elif curchoice == 2:
            self.upperPlot = "RDF"
        #event.Skip()
    def on_range(self, event):
        self.range = split('[,;:]?', self.text_ctrl_range.GetValue())
    def on_regular(self, event):
        self.text_ctrl_l0l1.Enable(self.radio_regular.GetValue())
    def on_l0l1(self, event):
        self.l0l1 = split('[,;:]?', self.text_ctrl_l0l1.GetValue())
    def on_backup(self, event):
        self.text_ctrl_l0l1.Enable(self.radio_regular.GetValue())
        file_choices = "Conf (*.txt)|*.txt"
        self.ConfPath = ""
        dlg = wx.FileDialog(self, message="Open configurational file...", defaultDir=getcwd(), defaultFile="", wildcard=file_choices, style=wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.ConfPath = dlg.GetPath()
            self.initConf = "2"
    def on_refreshSamp(self, event):
        try:
            f = int(self.text_ctrl_refreshSamp.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown Refresh Samples:", self.text_ctrl_refreshSamp.GetValue())
            return
        if f < 1:
            self.mainFrame.refresh_status_message("Too small value %d" % f + " and reset to 1")
            f = 1
            return
        self.lqfluid.RefreshSamp = f;    
        
    def on_stride(self, event):
        try:
            f = int(self.text_ctrl_stride.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown Stride:", self.text_ctrl_stride.GetValue())
            return
        if f < 1:
            self.mainFrame.refresh_status_message("Too small value %d" % f + " and reset to 1")
            f = 1
            return
        self.lqfluid.Stride = f;    
    def on_bau(self, event):
        try:
            self.bau = split('[,;]?', self.text_ctrl_bau.GetValue())
            beta = float(self.bau[0])
            alpha = float(self.bau[1])
            un = float(self.bau[2])
        except:
            self.mainFrame.refresh_status_message("Unknown beta alpha and U/N:" + +",".join(self.self.bau))
            return
        if beta < 0:
            self.mainFrame.refresh_status_message("Too small beta value %d" % beta)
            return
        self.lqfluid.dBeta = beta
        self.lqfluid.dAlpha = alpha
        self.lqfluid.dUn = un
    def on_lqPar(self, event):
        lqpar = split('[,;:]?', self.text_ctrl_lqPar.GetValue())
        try:
            self.lqfluid.eps = float(lqpar[0])
            self.lqfluid.sig = float(lqpar[1])
            self.lqfluid.rs = float(lqpar[2])
            self.lqfluid.rc = float(lqpar[3])
        except:
            self.mainFrame.refresh_status_message("Unknown Lennard-jones parameters:" + ",".join(lqpar))
            return
        
        self.lqfluid.ResetLQ(self.lqfluid.eps , self.lqfluid.sig , self.lqfluid.rs, self.lqfluid.rc)
        
    def record_figure(self):
        if self.record_figure_checked:
            figfile = "lqfluid%d_%d~%lg_%lg_%ld_%d_%ldplot%s.png" % (self.lqfluid.nDim0, self.lqfluid.nDim1, self.lqfluid.dDist, \
            self.lqfluid.dBeta, len(self.lqfluid.TR), self.lqfluid.Stride, self.lqfluid.uCurMove, strftime("%Y-%m-%d@%H~%M~%S"));
            self.fig.savefig(figfile) 
            
            print "save figure to:", figfile
        else:
            pass
