import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation


class Population:
    """
    A class to generate and evolve a lattice according to the rules of
    Greenberg-Hastings cellular automata for diffusion on excitable media.

    ...

    Attributes
    ----------
    act : int
        <desc>
    pas : int
        <desc>
    tau0 : int
        <desc>
    r : float
        <desc>
    size : int
        <desc>
    i_0 : float
        <desc>
    r_0 : float
        <desc>
    loc : float
        <desc>
    (s,i,r)count : float
        <desc>
    (s,i,r)time : list
        <desc>
    periodic : bool
        <desc>
    p : numpy.ndarray
        <desc>

    Methods
    -------
    count()
        <desc>
    nbr(i,j)
        <desc>
    check(loc)
        <desc>
    infect(loc)
        <desc>
        

    """
    def __init__(self,act=1,pas=1,r=1,size=1,i_0=0.0,r_0=None,periodic=False,p=None):
        """
        Parameters
        ----------
            act : int
                <desc>
            pas : int
                <desc>
            r : float
                <desc>
            size : int
                <desc>
            i_0 : float
                <desc>
            r_0 : float
                <desc>
            periodic : bool
                <desc>
            p : numpy.ndarray
                <desc>

        """
        self.act=act
        self.pas=pas
        self.tau0=self.act+self.pas

        self.loc = []
        self.scount,self.icount,self.rcount = 0,0,0
        self.stime,self.itime,self.rtime = [],[],[]

        if p is None:
            self.size = size
            self.i_0 = i_0

            if r_0 is None:
                self.r_0 = 0.5
            else:
                self.r_0 = r_0

##            self.p = np.random.randint(1,self.tau0+1,size=self.size*self.size)
            self.p = np.zeros(self.size*self.size,dtype=np.int8)
            endi0 = int(self.i_0*(self.size**2))
            endr0 = endi0+int((self.r_0*(1-self.i_0))*(self.size**2))

##            self.p[0:endi0] = 1
##            self.p[endi0:endr0] = self.act+1
            self.p[0:endi0] = np.random.randint(1,self.act+1)
            self.p[endi0:endr0] = np.random.randint(self.act+1,self.tau0+1)

            np.random.shuffle(self.p)
            self.p = self.p.reshape(self.size,self.size)
        else:
            self.p = p
            self.size = len(p[0])
            self.count()
            self.i_0 = self.icount/self.size**2
            self.r_0 = self.rcount/self.size**2

        self.periodic = periodic
        self.r = r
        self.default_nbh = [np.array([x,y]) for x in range(-self.size+1,self.size)
                            for y in range(-self.size+1,self.size)
                            if np.sqrt(x**2+y**2)<=r and [x,y]!=[0,0]]

    def count(self):
        """Counts the number of cells which are restive, active or passive.

        The attributes (s,i,r)count are the absolute numbers of the respective phases at the current time step,
        and the attributes (s,i,r)time are a list of the fractions of the respective phases for all time steps till the current one.

        Parameters
        ----------
        self : Population
            <desc>

        Returns
        ----------
        
        
        Raises
        ------

        """

        self.scount,self.icount,self.rcount = 0,0,0
        for i in range(self.size):
            for j in range(self.size):
                if self.p[i,j] == 0:
                    self.scount += 1
                elif self.p[i,j] in range(1,self.act+1):
                    self.icount += 1
                elif self.p[i,j] in range(self.act+1,self.act+self.pas+1):
                    self.rcount += 1
        self.stime.append(self.scount/self.size**2)
        self.itime.append(self.icount/self.size**2)
        self.rtime.append(self.rcount/self.size**2)
        return

    def nbr(self,i,j):
        """<desc>

        Parameters
        ----------
        <var> : <type>
            <desc>

        Returns
        ----------
        <var> : <type>
            <desc>
        
        Raises
        ------
        <err>
            <desc>
        """
        
        nbh = np.array([i,j])+self.default_nbh
        if self.periodic:
            c_nbrs = []
            for i in nbh%self.size:
                c_nbrs.append(self.p[tuple(i)])
        else:
            real_nbh = []
            for i in nbh:
                if 0<=i[0]<self.size and 0<=i[1]<self.size:
                    real_nbh.append(i)
            c_nbrs = []
            for i in real_nbh:
                c_nbrs.append(self.p[tuple(i)])
        return c_nbrs

    def check(self,loc):
        """<desc>

        Parameters
        ----------
        <var> : <type>
            <desc>

        Returns
        ----------
        <var> : <type>
            <desc>
        
        Raises
        ------
        <err>
            <desc>
        """

        self.loc = []
        nbhood = []
        for i in range(self.size):
            for j in range(self.size):
                if self.p[i,j] == 0:
                    if self.size==1:
                        continue
                    else:
                        nbhood = self.nbr(i,j)
                    for k in nbhood:
                        if k in range(1,self.act+1):
                            self.loc.append((i,j))
                            break
        return self.loc

    def infect(self,loc):
        """<desc>

        Parameters
        ----------
        <var> : <type>
            <desc>

        Returns
        ----------
        <var> : <type>
            <desc>
        
        Raises
        ------
        <err>
            <desc>
        """

        for i in range(self.size):
            for j in range(self.size):
                if self.p[i,j] >= self.tau0:
                    self.p[i,j] = 0
                if 1<=self.p[i,j]<self.tau0:
                    self.p[i,j] += 1
                if (i,j) in loc:
                    self.p[i,j] += 1
        return
####################################################################################

def run(com):
    """<desc>

    Parameters
    ----------
    <var> : <type>
        <desc>

    Returns
    ----------
    <var> : <type>
        <desc>
    
    Raises
    ------
    <err>
        <desc>
    """
    
##    for i in com:
##        i.count()
    for i in com:
        i.loc = i.check(i.loc)
    for i in com:
        i.infect(i.loc)
    return com

####################################################################################

def set_plot_text(pop,ax,cycle):
    core_len = int(len(pop))
    for txt in ax.texts:
        txt.set_text('');
    for i in range(core_len):
        for j in range(core_len):
            c='w' if 0<pop.T[i,j]<=cycle else 'k'
            s_txt = ax.text(i,j,'{0:2d}'.format(int(pop.T[i,j])),color=c,va='center',ha='center')
    return s_txt,

def plot(com,cid=None,ax=None,cbar=False,txt=False):
    core_len = int(len(com[0].p));cycle=int(com[0].tau0);n_states = cycle+1
    if ax == None:
        fig,ax = plt.subplots()
    ax.set_xticks(np.array(range(0,core_len))+0.5)
    ax.set_yticks(np.array(range(0,core_len))+0.5)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(axis='both', which='both', length=0)
    ax.set_aspect('equal')
    ax.grid(which='both')
    if cid:
        ax.set_title("$ID={}$".format(cid))
    if txt:
        s_txt, = set_plot_text(com[0].p,ax,cycle)

    cmap = colors.ListedColormap(['xkcd:pale grey','xkcd:darkish red','xkcd:almost black'])
    bounds = [0,0.5,com[0].act+0.5,com[0].tau0+0.5]
    norm = colors.BoundaryNorm(bounds,cmap.N)

    img = ax.imshow(com[0].p,cmap=cmap,norm=norm,interpolation='none')
    if cbar:
        if cbar == 'v':
            cbar = plt.colorbar(img,ax=ax,cmap=cmap,norm=norm,
                        boundaries=bounds,spacing='proportional',
                        ticks=range(n_states),orientation='vertical',
                        drawedges=False)
        elif cbar == 'h':
            cbar = plt.colorbar(img,ax=ax,cmap=cmap,norm=norm,
                        boundaries=bounds,spacing='proportional',
                        ticks=range(n_states),orientation='horizontal',
                        drawedges=False)
##        cbar.ax.tick_params(labelsize=20)
        fig.tight_layout()
        plt.show()
    return img,

def animate(com,data,txt=False,interval=100,T=50):
    core_len = int(len(com[0].p));cycle=int(com[0].tau0);n_states = cycle+1
    
    cmap = colors.ListedColormap(['xkcd:pale grey','xkcd:darkish red','xkcd:almost black'])
    bounds = [0,0.5,com[0].act+0.5,com[0].tau0+0.5]
    norm = colors.BoundaryNorm(bounds,cmap.N)        
        
    fig,ax = plt.subplots()
    img, = plot(com,ax=ax)

    if txt:
        txt, = set_plot_text(data[0],ax,cycle)
        def update(n,data):
            img.set_data(data[n])
            ax.set_title('$t={}$'.format(n))
            text, = set_plot_text(data[n],ax,cycle)
            return img,text,
    else:
        def update(n,data):
            img.set_data(data[n])
            ax.set_title('$t={}$'.format(n))
            return img,
    
    cax = fig.add_axes([0.83,0.11,0.03,0.775])
##        cax = fig.add_axes([0.05,0.05,0.9,0.03])
    cbar = plt.colorbar(img,cax=cax,cmap=cmap,norm=norm,
                        boundaries=bounds,spacing='proportional',
                        ticks=range(n_states),
                        drawedges=False,orientation='vertical')        

    ani = animation.FuncAnimation(fig, update, T, fargs=(data,),
                                  interval=interval, repeat=False)
    plt.show()
    return


if __name__=='__main__':
    pass
##    trans=0;obs=100;T=trans+obs
####    tau_i = 5 ; tau_r = 10 ; tau_0 = tau_i+tau_r
####    N = 1#tau_0
##
####   Population(self,i_0=None,r_0=None,p=None,size=None,
####              act=None,pas=None,r=1,periodic=False)
##
####    com=[Population(i_0=0.0,r_0=0.0)] ; com[0].p[5,5]=4;com[0].p[5,6]=9
##    com=[Population(size=4)] #; com[0].p[0,0]=1
##    N=com[0].size
##    plot(com,cbar='v',txt=True)
##    states=np.empty((T,int(N),int(N)),dtype=np.int8)
##    for t in range(T):
##        states[t] = np.copy(com[0].p)
##        com = run(com)
##    animate(com,states,T=T,txt=True)
