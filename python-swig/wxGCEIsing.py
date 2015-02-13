from mcss import Ising, PyJackknife_MeanStdErrTau_I

import wx
from numpy import zeros, int8, mean, std, arange, array
from re import split
from os import getcwd
from time import  time, strftime


import threading

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
import matplotlib.patches as patches
import matplotlib.transforms as transforms
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
import pylab

class CGCEIsing(Ising):
    def __init__(self):
        Ising.__init__(self)#must use 'self'
        self.RefreshSamp = 100
        self.Stride = 1
        self.Using_wolff = False
    def Reset(self, dim0, dim1):
        self.InitParameters()
        self.conf_t = zeros((dim0, dim1), int8)
        self.conf = self.PySetConf(self.conf_t)
        self.E = self.PyGetEnergyRange()#numpy array E(n) 
        self.H = zeros((len(self.E), 1), int) #for self.ising.PyTraj2Hist(ydata, self.ising.H)
        self.MH = zeros((self.nN + 1, 1), int) #for self.ising.PyTraj2Hist(ydata, self.ising.H)
        self.ETR = []#for energy
        self.MTR = []#for magnetism
        self.Stack_p, self.Seq_p = self.PySetGCEWolff()
    def ResetWolff(self):
        self.dGCEWF_AvgE = self.nCurHami
        self.dGCEWF_SegAccumE = 0
    def SetGCE(self, beta, alpha, un):
        self.dBeta = beta
        self.dAlpha = alpha
        self.dUn = un
    def RunGCE(self):
        if self.Using_wolff:
            for CurSample in range(self.RefreshSamp):
                self.WolffClusterGCE(self.Stride) #C++ extension Function from Ising in pymcsim model
                self.ETR.append(self.nCurHami)
                self.MTR.append(self.nCurMag)                
        else:
            Stride = self.Stride * self.nN
            for CurSample in range(self.RefreshSamp):
                self.MetroplisTrialGCE(Stride) #C++ extension Function from Ising in pymcsim model
                self.ETR.append(self.nCurHami)
                self.MTR.append(self.nCurMag)
                self.Gauge()
          
class FuncPanel(wx.Panel):
    """ A static panel for Ising
    """
    def __init__(self, parent, ID, mainFrame):
        self.ising = CGCEIsing()
        wx.Panel.__init__(self, parent, ID)
        self.mainFrame = mainFrame;
        self.label_l1l2 = wx.StaticText(self, -1, "L1,L2:  ", style=wx.ALIGN_CENTRE)
        self.text_ctrl_l1l2 = wx.TextCtrl(self, -1, "50,50")
        self.label_bau = wx.StaticText(self, -1, "Beta,Alpha,U/N:")
        self.text_ctrl_bau = wx.TextCtrl(self, -1, "0.431,1,-1.5")
        self.label_stride = wx.StaticText(self, -1, "Stride:")
        self.text_ctrl_stride = wx.TextCtrl(self, -1, "1")
        self.label_refreshSamp = wx.StaticText(self, -1, "Refresh Samples:")
        self.text_ctrl_refreshSamp = wx.TextCtrl(self, -1, "100")

        self.label_range = wx.StaticText(self, -1, "From:To ", style=wx.ALIGN_CENTRE)
        self.text_ctrl_range = wx.TextCtrl(self, -1, "-2000:")
        
        self.checkbox_wolff = wx.CheckBox(self, -1, "Wolff Cluster")
        self.label_wolff = wx.StaticText(self, -1, "Average Length:")
        self.text_ctrl_wolff_Eavg = wx.TextCtrl(self, -1, "10")
        
        self.radio_box_initial = wx.RadioBox(self, -1, "Initial state", choices=["Ground", "Random", "Anti-ground", "Backup file"], majorDimension=4, style=wx.RA_SPECIFY_ROWS)
        self.radio_box_leftplot = wx.RadioBox(self, -1, "Left subplot", choices=["Energy Trajectory", "Energy Histogram", "Magnetism Trajectory", "Magnetism Histogram"], majorDimension=4, style=wx.RA_SPECIFY_ROWS)
        self.button_reset = wx.Button(self, -1, "Reset")
        self.button_start = wx.Button(self, -1, "Start")
        self.button_analysis = wx.Button(self, -1, "Analysis")
        self.__do_layout()

        self.build_frame(parent)
        
        self.reset = False;
        self.ploting = False
        self.keeping = False
        self.job = None     
        self.record_figure_checked = False   
        self.transX = matplotlib.transforms.blended_transform_factory(self.axes_left.transData, self.axes_left.transAxes)    
        self.transY = matplotlib.transforms.blended_transform_factory(self.axes_left.transAxes, self.axes_left.transData)
        
        self.on_initial(None)
        self.on_leftplot(None)
        #self.cmap = matplotlib.cm.colors.ListedColormap(['b', 'y'])
        self.redraw_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_redraw_timer, self.redraw_timer)        
        self.redraw_timer.Start(400)

    def __do_layout(self):
        grp1_box = wx.StaticBox(self, -1, "Ising model")
        grp1_sizer = wx.StaticBoxSizer(grp1_box, wx.HORIZONTAL)
        grp1_box.SetMinSize((220, 70))
        
        sizer_model1 = wx.BoxSizer(wx.VERTICAL)
        sizer_model1.Add(self.label_l1l2, 0, wx.ALL, 2)
        sizer_model1.Add(self.text_ctrl_l1l2, 0, wx.ALL, 2)
        sizer_model1.Add(self.label_bau, 0, wx.ALL, 2)
        sizer_model1.Add(self.text_ctrl_bau, 0, wx.ALL, 2)

        sizer_model2 = wx.BoxSizer(wx.VERTICAL)
        sizer_model2.Add(self.label_stride, 0, wx.ALL, 2)
        sizer_model2.Add(self.text_ctrl_stride, 0, wx.ALL, 2)
        self.text_ctrl_stride.Bind(wx.EVT_TEXT, self.on_stride)
        sizer_model2.Add(self.label_refreshSamp, 0, wx.ALL, 2)
        sizer_model2.Add(self.text_ctrl_refreshSamp, 0, wx.ALL, 2)
        self.text_ctrl_refreshSamp.Bind(wx.EVT_TEXT, self.on_refreshSamp)
               
        grp1_sizer.Add(sizer_model1, 0, wx.ALL, 2)
        grp1_sizer.Add(sizer_model2, 0, wx.ALL, 2)

        
        grp2_sizer = wx.BoxSizer(wx.VERTICAL)
        grp2_sizer.Add(self.radio_box_initial, 0, wx.ALL, 1)
        self.radio_box_initial.Bind(wx.EVT_RADIOBOX, self.on_initial)
        self.radio_box_initial.SetSelection(1)
        self.radio_box_initial.SetMinSize((116, 130))

        grp3_sizer = wx.BoxSizer(wx.VERTICAL)
        grp3_sizer.Add(self.radio_box_leftplot, 0, wx.ALL, 1)
        self.radio_box_leftplot.Bind(wx.EVT_RADIOBOX, self.on_leftplot)
        self.radio_box_leftplot.SetSelection(0)
        self.radio_box_leftplot.SetMinSize((175, 130))
            
        grp5_box = wx.StaticBox(self, -1, "")
        grp5_sizer = wx.StaticBoxSizer(grp5_box, wx.VERTICAL)
        grp5_box.SetMinSize((158, 60))
        grp5_sizer.Add(self.checkbox_wolff, 0, wx.ALL, 0)
        grp5_sizer.Add(self.label_wolff, 0, wx.ALL, 1)
        grp5_sizer.Add(self.text_ctrl_wolff_Eavg, 0, wx.ALL, 0)
        self.text_ctrl_wolff_Eavg.Bind(wx.EVT_TEXT, self.on_wolff_Eavg)
        self.text_ctrl_wolff_Eavg.Disable()
        self.text_ctrl_wolff_Eavg.SetMinSize((80, 24))
        self.checkbox_wolff.Bind(wx.EVT_CHECKBOX, self.on_check_wolff)
        
        grpTraj_box = wx.StaticBox(self, -1, "Trajectory Range")
        grpTraj_sizer = wx.StaticBoxSizer(grpTraj_box, wx.HORIZONTAL)
        grpTraj_box.SetMinSize((158, 38))      
        grpTraj_sizer.Add(self.label_range, 0, wx.ALL, 4)
        grpTraj_sizer.Add(self.text_ctrl_range, 0, wx.ALL, 1)
        self.text_ctrl_range.SetMinSize((80, 25))
        self.text_ctrl_range.Bind(wx.EVT_TEXT, self.on_range)
        
        grpWJ_sizer = wx.BoxSizer(wx.VERTICAL)
        grpWJ_sizer.Add(grp5_sizer, 0, wx.ALL, 1)
        grpWJ_sizer.Add(grpTraj_sizer, 0, wx.ALL, 1)
 
        grp6_box = wx.StaticBox(self, -1, "")
        grp6_sizer = wx.StaticBoxSizer(grp6_box, wx.VERTICAL)
        grp6_box.SetMinSize((100, 130))
        grp6_sizer.Add(self.button_reset, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        grp6_sizer.Add(self.button_start, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        grp6_sizer.Add(self.button_analysis, border=4, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.button_reset.Bind(wx.EVT_BUTTON, self.on_reset)
        self.button_start.Bind(wx.EVT_BUTTON, self.on_start)
        self.button_analysis.Bind(wx.EVT_BUTTON, self.on_analysis)
        self.button_start.Disable()  
        self.button_analysis.Disable()
        self.sizer_body = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer_body.Add(grp1_sizer, border=5, flag=wx.ALL)
        self.sizer_body.Add(grp2_sizer, border=5, flag=wx.ALL)
        self.sizer_body.Add(grp3_sizer, border=5, flag=wx.ALL)
        #self.sizer_body.AddSpacer(24)
        self.sizer_body.Add(grpWJ_sizer, border=5, flag=wx.ALL)
        self.sizer_body.Add(grp6_sizer, border=5, flag=wx.ALL)
        self.SetSizer(self.sizer_body)
        self.sizer_body.Fit(self)
    def build_frame(self, mainPanel):
        self.fig = Figure((8.0, 3.5), dpi=100)
        self.axes_left = self.fig.add_subplot(121)
        self.axes_right = self.fig.add_subplot(122)
        self.axes_left.set_axis_bgcolor('gray')
        self.axes_right.set_axis_bgcolor('gray')
        self.fig.subplots_adjust(bottom=0.14, left=0.11, right=0.96)
        self.canvas = FigCanvas(mainPanel, -1, self.fig)
        # Create the navigation toolbar, tied to the canvas
        self.toolbar = NavigationToolbar(self.canvas)
        # Layout with box sizers
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        #self.vbox.AddSpacer(10)        
        self.vbox.Add(self, 0, flag=wx.ALL | wx.EXPAND)
        mainPanel.SetSizer(self.vbox)
        self.vbox.Fit(self.mainFrame)

    """Event functions"""
    def on_reset(self, event):
        self.reset = False;
        self.button_start.Disable()
        l1l2 = split('[,;]?', self.text_ctrl_l1l2.GetValue())
        bau = split('[,;]?', self.text_ctrl_bau.GetValue())
        self.range = split('[,;:]?', self.text_ctrl_range.GetValue())
        #print l1l2q, bau, self.range 
        if len(l1l2) != 2:
            self.mainFrame.refresh_status_message("Unknown input:" + self.text_ctrl_l1l2.GetValue())
            return False
        if len(bau) != 3:
            self.mainFrame.refresh_status_message("Unknown input:" + self.self.text_ctrl_bau.GetValue())
            return False
        self.ising.Reset(int(l1l2[0]), int(l1l2[1]))
        initcf = int(self.initConf)
        if initcf == 2:
            if not self.ising.LoadConf(self.ConfPath):
                self.mainFrame.refresh_status_message("Load conf file error:" + self.ConfPath)     
                return False     
        else:
            self.ising.InitConf(initcf)
        self.ising.Gauge() 
        self.ising.SetGCE(float(bau[0]), float(bau[1]), float(bau[2]))
        
        self.on_check_wolff(None)
        self.on_refreshSamp(None)
        self.on_stride(None)

        
        self.button_start.Enable()
                
        self.axes_left.clear()
        self.axes_left.text(0.1, 0.85, "Generalized canonical ensemble", color="blue", size=12, transform=self.axes_left.transAxes)
        self.axes_left.text(0.1, 0.6, "$Z=\sum_E g(E) e ^{-\\beta E-\\alpha(E-U)^2/(2N)}$", color="red", size=14, transform=self.axes_left.transAxes)
        self.axes_left.text(0.1, 0.4, "Ising model $H=-\sum_{<i,j>} s_i \cdot s_j$", color="red", size=14, transform=self.axes_left.transAxes)
        self.axes_left.set_axis_bgcolor('black')
        # for configuration plotting
        self.axes_right.clear()
        self.axes_right.set_title(r'Ising Model $N=%d\times%d$, $bright=\downarrow$, $green=\uparrow$' % (self.ising.nDim0, self.ising.nDim1), size=10)
        self.axes_right.set_xlabel('L2', fontsize=9)
        self.axes_right.set_ylabel('L1', fontsize=9)
        self.axes_right.grid(False)# interpolation='nearest', animated=True,matplotlib.cm.YlGn
        self.im = self.axes_right.imshow(self.ising.conf, cmap=matplotlib.cm.BuGn, interpolation='nearest', origin='lower', animated=True)

        pylab.setp(self.axes_right.get_xticklabels(), fontsize=9)
        pylab.setp(self.axes_right.get_yticklabels(), fontsize=9)
        """for tick in self.axes_right.yaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = True
            #tick.label2.set_color('green')
        """
        self.canvas.draw()  
        
        self.reset = True;
        return True
                 
    def on_start(self, event):
        if not self.keeping:
            self.keeping = True;
            label = "Pause"
            self.button_reset.Disable()
            self.checkbox_wolff.Disable()
            self.button_analysis.Disable()
            self.text_ctrl_l1l2.Disable()
            self.text_ctrl_bau.Disable()
            self.radio_box_initial.Disable()
            self.starting_time = time()#clock() takes negative when it lasts long  
        else:
            self.keeping = False;
            label = "Start"
            self.button_reset.Enable()
            self.checkbox_wolff.Enable()
            self.button_analysis.Enable()
            self.text_ctrl_l1l2.Enable()
            self.text_ctrl_bau.Enable()
            self.radio_box_initial.Enable()

        self.button_start.SetLabel(label)
        
    
    def on_redraw_timer(self, event):
        if self.ploting:
            return
        if not self.keeping:
            return
        
        if  self.job == None:
            self.job = threading.Thread(target=self.ising.RunGCE)# args=(self.ising, self.u)
            self.job.start()
        else:
            if not self.job.isAlive():
                self.ploting = True
                self.draw_figure()
                self.mainFrame.refresh_status_message("%lg seconds" % (time() - self.starting_time))
                self.job = threading.Thread(target=self.ising.RunGCE)
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
    def draw_figure(self):
        """ Redraws the figure
        """
        self.axes_left.clear()
        #self.axes_left.lines.remove(self.axes_left_line)
        self.axes_left.set_title(r'GCE %s: $\beta(E)=\beta+\alpha(E-U)/N$' % ("Wolff" if self.ising.Using_wolff else ""), size=10)
        self.axes_left.grid(True, color='gray')
        self.axes_left.set_axis_bgcolor('black')

        if self.LeftPlot == 'ETraj':
            if not self.get_data(self.ising.ETR):
                return
            self.axes_left.set_xlabel('Samples', fontsize=9)
            self.axes_left.set_ylabel('Energy', fontsize=9)
            self.axes_left.plot(self.tdata, self.ydata)

        elif self.LeftPlot == 'EHist':
            if not self.get_data(self.ising.ETR):
                return
            self.axes_left.set_xlabel('Energy', fontsize=9)
            self.axes_left.set_ylabel('Histogram', fontsize=9)
            e0, e1, Erange = self.ising.PyETraj2Hist(self.ydata, self.ising.H)
            self.axes_left.plot(arange(e0 / 4, e1 / 4 + 1) * 4, Erange[(e0 + 2 * self.ising.nN) / 4:(e1 + 2 * self.ising.nN) / 4 + 1], '.', markersize=3)
        elif self.LeftPlot == 'MTraj':
            if not self.get_data(self.ising.MTR):
                return
            self.axes_left.set_xlabel('Samples', fontsize=9)
            self.axes_left.set_ylabel('Magnetism', fontsize=9)
            self.axes_left.plot(self.tdata, self.ydata, 'g-')
        elif self.LeftPlot == 'MHist':
            if not self.get_data(self.ising.MTR):
                return
            self.axes_left.set_xlabel('Magnetism', fontsize=9)
            self.axes_left.set_ylabel('Histogram', fontsize=9)
            m0, m1, Hrange = self.ising.PyMTraj2Hist(self.ydata, self.ising.MH)
            m0i = (m0 + self.ising.nN) / 2
            m1i = (m1 + self.ising.nN) / 2
            self.axes_left.plot(arange(m0 , m1 + 2, 2) , Hrange[m0i:m1i + 1 ], 'g.', markersize=2)

        for label in self.axes_left.xaxis.get_ticklabels():
            # label is a Text instance
            #label.set_color('red')
            label.set_rotation(20)
            label.set_fontsize(9)
        #pylab.setp(self.axes_left.get_xticklabels(), fontsize=9)
        pylab.setp(self.axes_left.get_yticklabels(), fontsize=9)

        #plotting configuration  
        
        if self.ising.nInitconf == -1:#YlGn   
            self.axes_right.imshow(self.ising.conf, cmap=matplotlib.cm.BuGn, interpolation='nearest', origin='lower', animated=True)
        else:
            self.im.set_array(self.ising.conf) #not suitable for ground state
 
        self.canvas.draw()  
        self.record_figure()
    def on_analysis(self, event):
        if  len(self.ydata) == 0:
            return

        if self.LeftPlot == 'ETraj':
            self.draw_ananlysis(r'$\bar E=%g\pm%g$, $\tau_i=%g$', Y=True)
        elif self.LeftPlot == 'EHist':
            self.draw_ananlysis(r'$\bar E=%g\pm%g$, $\tau_i=%g$', Y=False)
        elif self.LeftPlot == 'MTraj':
            self.draw_ananlysis(r'$\bar M=%g\pm%g$, $\tau_i=%g$', Y=True)
        elif self.LeftPlot == 'MHist':
            self.draw_ananlysis(r'$\bar M=%g\pm%g$, $\tau_i=%g$', Y=False)
        return
    def draw_ananlysis(self, labels, Y):
            m, s, merr , ptau, taui, tauCorr = PyJackknife_MeanStdErrTau_I(self.ydata, 10, 80)
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
  
            self.axes_left.text(0.05, 0.9, labels % (m, merr, ptau), horizontalalignment='left', verticalalignment='top', color="red", size=10, transform=self.axes_left.transAxes)
            self.axes_left.lines.extend([l1])
            self.axes_left.add_patch(rect)
            self.canvas.draw() 
            self.record_figure()
    def on_initial(self, event):
        self.curchoice = self.radio_box_initial.GetSelection()
        if self.curchoice == 0:
            self.initConf = "-1"
        elif self.curchoice == 1:
            self.initConf = "0"
        elif self.curchoice == 2:
            self.initConf = "1"
        else:
            file_choices = "Conf (*.txt)|*.txt"
            self.ConfPath = ""
            dlg = wx.FileDialog(
                self,
                message="Open configurational file...",
                defaultDir=getcwd(),
                defaultFile="",
                wildcard=file_choices,
                style=wx.OPEN)
        
            if dlg.ShowModal() == wx.ID_OK:
                self.ConfPath = dlg.GetPath()
                self.initConf = "2"
        #event.Skip()
  
    def on_leftplot(self, event):
        curchoice = self.radio_box_leftplot.GetSelection()
        if curchoice == 0:
            self.LeftPlot = "ETraj"
        elif curchoice == 1:
            self.LeftPlot = "EHist"
        elif curchoice == 2:
            self.LeftPlot = "MTraj"
        elif curchoice == 3:
            self.LeftPlot = "MHist"
        #event.Skip()
    def on_range(self, event):
        self.range = split('[,;:]?', self.text_ctrl_range.GetValue())
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
        self.ising.RefreshSamp = f;    
        
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
        self.ising.Stride = f;    
    def on_check_wolff(self, event):
        self.text_ctrl_wolff_Eavg.Enable(self.checkbox_wolff.IsChecked())
        self.ising.Using_wolff = self.checkbox_wolff.IsChecked()
        if self.ising.Using_wolff:
            self.on_wolff_Eavg(None)
            self.ising.ResetWolff()
    def on_wolff_Eavg(self, event):
        try:
            f = int(self.text_ctrl_wolff_Eavg.GetValue())
        except:
            self.mainFrame.refresh_status_message("Unknown wolff_Eavg:", self.text_ctrl_wolff_Eavg.GetValue())
            return
        if f < 2:
            self.mainFrame.refresh_status_message("Too small value %d" % f + " and reset to 2")
            f = 2
            return
        self.ising.nGCEWF_AvgC = f;      
        
    """for Mainframe calling functions"""
    def message_dlg(self, titile, msg):
            dlg = wx.MessageDialog(self, msg, titile, wx.OK | wx.ICON_INFORMATION)
            dlg.ShowModal() 
            dlg.Destroy()
            
    def save_traj(self):
        if self.keeping: 
            self.message_dlg("Error", "Stop runing first!")
            return
        if not self.reset: 
            self.message_dlg("Error", "No resetting!")
            return
        
        file_choices = "traj (*.txt)|*.txt"
        if self.ising.Using_wolff:
            trajfile = "gcewfp%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ldTraj.txt" % (self.ising.nDim0, self.ising.nDim1, ord(self.ising.cQ), self.ising.nInitconf, \
            self.ising.dBeta, self.ising.dAlpha, self.ising.dUn, \
            self.ising.Stride, self.ising.RefreshSamp);
        else:
            trajfile = "gce%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ldTraj.txt" % (self.ising.nDim0, self.ising.nDim1, ord(self.ising.cQ), self.ising.nInitconf, \
            self.ising.dBeta, self.ising.dAlpha, self.ising.dUn, \
            self.ising.Stride, self.ising.RefreshSamp);
                        
        dlg = wx.FileDialog(self, message="Save trajectory file...",
                defaultDir=getcwd(), defaultFile=trajfile, wildcard=file_choices,
                style=wx.SAVE)
                
        if dlg.ShowModal() == wx.ID_OK:
                trajfile = dlg.GetPath()
                f = open(trajfile, 'w')
                if f:
                    for item in self.ising.ETR:
                        f.write(str(item) + "\n")

    def save_conf(self):
        if self.keeping: 
            self.message_dlg("Error", "Stop runing first!")
            return
        if not self.reset: 
            self.message_dlg("Error", "No resetting!")
            return 
        
        file_choices = "Conf (*.conf)|*.conf"
        if self.ising.Using_wolff:
            trajfile = "gcewfi%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ld.conf" % (self.ising.nDim0, self.ising.nDim1, self.ising.nInitconf, \
            self.ising.dBeta, self.ising.dAlpha, self.ising.dUn, \
            self.ising.Stride, self.ising.RefreshSamp);
        else:
            trajfile = "gcei%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ld.conf" % (self.ising.nDim0, self.ising.nDim1, self.ising.nInitconf, \
            self.ising.dBeta, self.ising.dAlpha, self.ising.dUn, \
            self.ising.Stride, self.ising.RefreshSamp);
                        
        dlg = wx.FileDialog(self, message="Save configuration file...",
                defaultDir=getcwd(), defaultFile=trajfile, wildcard=file_choices,
                style=wx.SAVE)
                
        if dlg.ShowModal() == wx.ID_OK:
                trajfile = dlg.GetPath()
                print trajfile
                self.ising.DumpConf(trajfile)

    def record_figure(self):
        if self.record_figure_checked:
            figfile = "gcei%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ld_%ldplot%s.png" % (self.ising.nDim0, self.ising.nDim1, self.ising.nInitconf, \
            self.ising.dBeta, self.ising.dAlpha, self.ising.dUn, \
            self.ising.Stride, self.ising.RefreshSamp, self.ising.uCurMove, strftime("%Y-%m-%d@%H~%M~%S"));
            self.fig.savefig(figfile) 
            print "save figure to:", figfile
        else:
            pass
