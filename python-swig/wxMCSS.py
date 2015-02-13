import wx
import wxGCEIsing
import wxWLXIsing
import wxGCEPotts
import wxWLXPotts
import wxMulticaPotts
import wxSingleParticle
import wxOffLat2D
import wxOffLat3D


class MainFrame(wx.Frame):
    """ The main frame of the application
    """
    title = 'Statistical model simulator'
    
    def __init__(self):
        wx.Frame.__init__(self, None, -1, self.title)
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas 
             * mpl navigation toolbar
             * Model panel
             * status bar
        """
        self.panel = wx.Panel(self)
        self.funcpanl = wxWLXIsing.FuncPanel(self.panel, -1, self);#BoundControlBox(self.panel, -1, "X min", 0)

    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        self.m_save_traj = menu_file.Append(-1, "&Save trajectory\tCtrl-T", "Save trajectory to file")
        self.Bind(wx.EVT_MENU, self.on_save_traj, self.m_save_traj)
        self.m_save_conf = menu_file.Append(-1, "&Save Configuration\tCtrl-F", "Save configuration to file")
        self.Bind(wx.EVT_MENU, self.on_save_conf, self.m_save_conf)
        self.m_record_figure = menu_file.Append(-1, "&Record Figures\tCtrl-R", "Record plotting figures to file", wx.ITEM_CHECK)
        self.Bind(wx.EVT_MENU, self.on_record_figure, self.m_record_figure)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
   
        menu_model = wx.Menu()
       
        m_ising = wx.Menu()
        self.m_ising_gce_wf = wx.MenuItem(m_ising, wx.NewId(), "&Generalized Canonical Ensemble\tF2", "", wx.ITEM_NORMAL)
        m_ising.AppendItem(self.m_ising_gce_wf)
        self.m_ising_wl = wx.MenuItem(m_ising, wx.NewId(), "&Wang-Landua Method\tF3", "", wx.ITEM_NORMAL)
        m_ising.AppendItem(self.m_ising_wl)
        self.m_ising_cell_auto = wx.MenuItem(m_ising, wx.NewId(), "&Cellular Automata", "Cellular automata on Ising model", wx.ITEM_NORMAL)
        m_ising.AppendItem(self.m_ising_cell_auto)
        menu_model.AppendMenu(wx.NewId(), "Ising", m_ising, "Using 2D Ising model")
        
        m_potts = wx.Menu()
        self.m_potts_gce_wf = wx.MenuItem(m_potts, wx.NewId(), "&Generalized Canonical Ensemble\tF4", "", wx.ITEM_NORMAL)
        m_potts.AppendItem(self.m_potts_gce_wf)
        self.m_potts_wl = wx.MenuItem(m_potts, wx.NewId(), "&Wang-Landau Method\tF5", "", wx.ITEM_NORMAL)
        m_potts.AppendItem(self.m_potts_wl)
        self.m_potts_multican = wx.MenuItem(m_potts, wx.NewId(), "&Multicanonical Method\tF6", "", wx.ITEM_NORMAL)
        m_potts.AppendItem(self.m_potts_multican)
        self.m_potts_ge = wx.MenuItem(m_potts, wx.NewId(), "&Generalized Ensemble Method", "", wx.ITEM_NORMAL)
        m_potts.AppendItem(self.m_potts_ge)        
        menu_model.AppendMenu(wx.NewId(), "Potts", m_potts, "Using 2D Q-state Potts model")
                
        m_fluid = wx.Menu()
        self.m_fld_latoff2_mc = wx.MenuItem(m_fluid, wx.NewId(), "Off-Lattice 2D Model MC\tF7", "", wx.ITEM_NORMAL)
        m_fluid.AppendItem(self.m_fld_latoff2_mc)
        self.m_fld_latoff2_md = wx.MenuItem(m_fluid, wx.NewId(), "Off-Lattice 2D Model MD", "", wx.ITEM_NORMAL)
        m_fluid.AppendItem(self.m_fld_latoff2_md)
        self.m_fld_latoff3_mc = wx.MenuItem(m_fluid, wx.NewId(), "Off-Lattice 3D Model MC", "", wx.ITEM_NORMAL)
        m_fluid.AppendItem(self.m_fld_latoff3_mc)
        self.m_fld_latoff3_md = wx.MenuItem(m_fluid, wx.NewId(), "Off-Lattice 3D Model MD", "", wx.ITEM_NORMAL)
        m_fluid.AppendItem(self.m_fld_latoff3_md)
        menu_model.AppendMenu(wx.NewId(), "Fluid", m_fluid, "")
    
        m_polymer = wx.Menu()
        self.m_ply_hp_mc = wx.MenuItem(m_polymer, wx.NewId(), "HP Model MC", "", wx.ITEM_NORMAL)
        m_polymer.AppendItem(self.m_ply_hp_mc)
        menu_model.AppendMenu(wx.NewId(), "Polymer", m_polymer, "protein folding")
            
        menu_model.AppendSeparator()
        m_sig_pfun = menu_model.Append(-1, "Single Particle\tF8", "1D Potential function")
        self.Bind(wx.EVT_MENU, self.on_menu_sig_pfun, m_sig_pfun)
        m_selfavoiding = menu_model.Append(-1, "Self avoiding walker", "4D Snake")
        self.Bind(wx.EVT_MENU, self.on_selfavoiding, m_selfavoiding)
        
        self.Bind(wx.EVT_MENU, self.on_menu_ising_gce_wf, self.m_ising_gce_wf)
        self.Bind(wx.EVT_MENU, self.on_menu_ising_wl, self.m_ising_wl)
        self.Bind(wx.EVT_MENU, self.on_menu_ising_cell_auto, self.m_ising_cell_auto)
        self.Bind(wx.EVT_MENU, self.on_menu_potts_gce_wf, self.m_potts_gce_wf)
        self.Bind(wx.EVT_MENU, self.on_menu_potts_wl, self.m_potts_wl)
        self.Bind(wx.EVT_MENU, self.on_menu_potts_multican, self.m_potts_multican)
        self.Bind(wx.EVT_MENU, self.on_menu_fld_latoff2_mc, self.m_fld_latoff2_mc)
        self.Bind(wx.EVT_MENU, self.on_menu_fld_latoff3_mc, self.m_fld_latoff3_mc)

        # end wxGlade

        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)      

        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_model, "&Model")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)
        
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def on_save_traj(self, event):
        if hasattr(self.funcpanl, 'save_traj'):
            self.funcpanl.on_save_traj()
        else:
            pass
    def on_save_conf(self, event):
        if hasattr(self.funcpanl, 'save_conf'):
            self.funcpanl.on_save_conf()
        else:
            pass
    def on_record_figure(self, event):
        self.funcpanl.record_figure_checked = self.m_record_figure.IsChecked()
        print "Record figure:", self.m_record_figure.IsChecked()
    def on_menu_ising_gce_wf(self, event): 
        self.panel.Destroy()#delete self.funcpanl automatically
        self.panel = wx.Panel(self)
        self.funcpanl = wxGCEIsing.FuncPanel(self.panel, -1, self); 
        self.panel.Fit()
        self.statusbar.SetStatusText('')
    def on_menu_ising_wl(self, event):
        self.panel.Destroy()#delete self.funcpanl automatically
        self.panel = wx.Panel(self)
        self.funcpanl = wxWLXIsing.FuncPanel(self.panel, -1, self); 
        self.panel.Fit()
        self.statusbar.SetStatusText('')
    def on_menu_ising_cell_auto(self, event): 
        print "Event handler `on_menu_ising_cell_auto' not implemented!"

    def on_menu_potts_gce_wf(self, event): 
        self.panel.Destroy()#delete self.funcpanl automatically
        self.panel = wx.Panel(self)
        self.funcpanl = wxGCEPotts.FuncPanel(self.panel, -1, self); 
        self.panel.Fit()
        self.statusbar.SetStatusText('')

    def on_menu_potts_wl(self, event): 
        self.panel.Destroy()#delete self.funcpanl automatically
        self.panel = wx.Panel(self)
        self.funcpanl = wxWLXPotts.FuncPanel(self.panel, -1, self); 
        self.panel.Fit()
        self.statusbar.SetStatusText('')

    def on_menu_potts_multican(self, event): 
        self.panel.Destroy()#delete self.funcpanl automatically
        self.panel = wx.Panel(self)
        self.funcpanl = wxMulticaPotts.FuncPanel(self.panel, -1, self); 
        self.panel.Fit()
        self.statusbar.SetStatusText('')
    def on_menu_sig_pfun(self, event):  
        self.panel.Destroy()#delete self.funcpanl automatically
        self.panel = wx.Panel(self)
        self.funcpanl = wxSingleParticle.FuncPanel(self.panel, -1, self); 
        self.panel.Fit()
        self.statusbar.SetStatusText('')

    def on_selfavoiding(self, event):  
        pass
    def on_menu_fld_latoff2_mc(self, event):
        self.panel.Destroy()#delete self.funcpanl automatically
        self.panel = wx.Panel(self)
        self.funcpanl = wxOffLat2D.FuncPanel(self.panel, -1, self); 
        self.panel.Fit()
        self.statusbar.SetStatusText('')

    def on_menu_fld_latoff3_mc(self, event):
        self.panel.Destroy()#delete self.funcpanl automatically
        self.panel = wx.Panel(self)
        self.funcpanl = wxOffLat3D.FuncPanel(self.panel, -1, self); 
        self.panel.Fit()
        self.statusbar.SetStatusText('')
                  
    def on_exit(self, event):
        self.Destroy()
        
    def on_about(self, event):
        msg = """ Statistical model simulator based on wxPython with matplotlib:
         * 2009, Dec 19
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
    
    def refresh_status_message(self, msg, refresh_len_ms=1500):
        self.statusbar.SetStatusText(msg)
        """
        self.timeroff = wx.Timer(self)
        self.Bind(
            wx.EVT_TIMER,
            self.on_refresh_status_off,
            self.timeroff)
        self.timeroff.Start(refresh_len_ms, oneShot=True)
        """    
    def on_refresh_status_off(self, event):
        self.statusbar.SetStatusText('')

if __name__ == '__main__':
    app = wx.PySimpleApp()
    app.frame = MainFrame()
    app.frame.Show()
    app.MainLoop()


