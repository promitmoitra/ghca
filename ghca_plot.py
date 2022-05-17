import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import ghca_main as ca
import ghca_core as cacore


def run_config_id(pair,core,config_id,T=1,emb=True,return_str=False):
    """

    """
    active = pair[0] ; passive = pair[1] ; n_states = active+passive+1
    str_state = cacore.base_conv(config_id,n_states);str_state = str_state.rjust(core,'0')
    
    if emb:
        grid_size = np.sqrt(core)+4
        p = np.zeros((int(grid_size),int(grid_size)))
        q = cacore.str_to_state(str_state,n_states)
        if int(np.sqrt(core))%2==0:
            p[int(grid_size/2)-int(np.sqrt(core)/2):int(grid_size/2)+int(np.sqrt(core)/2),
            int(grid_size/2)-int(np.sqrt(core)/2):int(grid_size/2)+int(np.sqrt(core)/2)]=q
        else:
            p[int(grid_size/2)-int(np.sqrt(core)/2):int(grid_size/2)+int(np.sqrt(core)/2)+1,
            int(grid_size/2)-int(np.sqrt(core)/2):int(grid_size/2)+int(np.sqrt(core)/2)+1]=q
    else:
        grid_size = np.sqrt(core)        
        p = cacore.str_to_state(str_state,n_states)

    com = [ca.Population(p=np.copy(p),act=active,pas=passive)]

    states=np.empty((T,int(grid_size),int(grid_size)),dtype=np.int8);str_states = []
    for t in range(T):
        states[t] = np.copy(np.array(com[0].p,dtype=np.int8))        
        state = cacore.state_to_str(com[0].p,n_states)
        str_states.append(state)
        com = ca.run(com)

    if return_str:
        return str_states,p
    else:
        return states,p


def plot_perst_prob(core,flag,anchor=1,show=True):
    """
    
    """
    cycle_set = list(range(1,18))
    fixed_len_pairs = [(i,j) for i in cycle_set for j in cycle_set if i+j==cycle]
    fixed_passive = anchor ; inc_act_pairs = [(i,fixed_passive) for i in cycle_set]
    fixed_active = anchor ; inc_pas_pairs = [(fixed_active,i) for i in cycle_set]

    plt.ylim(-0.1,1.1)
    plt.ylabel("$P_{perst}$")

    if flag=='Fixed\ cycle':
        param_pair=fixed_len_pairs
        plt.xticks(range(len(param_pair)),labels=param_pair)
        plt.xlabel("$States\ -\ (Active,Passive)$")
        plt.title("$Core\ Size:\ {0},\ Cycle\ Length:\ {1}$".format(core,cycle))
    elif flag=='Fixed\ passive':
        param_pair=inc_act_pairs
        plt.xticks(range(len(param_pair)),labels=[i[0] for i in param_pair])
        plt.xlabel("$Active$")
        plt.title("$Core\ Size:\ {0}$".format(core))
    elif flag=='Fixed\ active':
        param_pair=inc_pas_pairs
        plt.xticks(range(len(param_pair)),labels=[i[1] for i in param_pair])
        plt.xlabel("$Passive$")
        plt.title("$Core\ Size:\ {0}$".format(core))
    else:
        print("Incorrect flag!")

    act_frac = {}
    for pair in param_pair:
        act_config_ids = np.load(data_path+"act_config_ids_states-({0:02d},{1:02d})_core-{2:02d}.npy".format(*pair,core))
        frac = len(act_config_ids)/(sum(pair)+1)**core
        act_frac[pair]=frac

    act_conf_frac = list(act_frac.values())
    plt.plot(range(len(param_pair)),act_conf_frac,'-o',label="${0:01d}$".format(anchor))
    plt.legend(loc='best',title="${}$".format(flag))
    if show:
        plt.show()
    return


def plot_param_space(core):
    """
    Fraction of init configs that persist (for given core size) = basin[passive-1,active-1]

    """
    max_cyc = 18 
    len_set = list(range(1,max_cyc))
    param_pairs=[(i,j) for i in len_set for j in len_set]
    param_space = np.zeros((max_cyc-1,max_cyc-1))
    for pair in param_pairs:
        act = pair[0];pas=pair[1]
        conf_id = np.load(data_path+"act_config_ids_states-({0:02d},{1:02d})_core-{2:02d}.npy".format(act,pas,core))
        param_space[pas-1,act-1] = len(conf_id)/(np.sum(pair)+1)**core
    plt.xticks(range(max_cyc-1),labels=len_set);plt.yticks(range(max_cyc-1),labels=len_set)
    plt.xlabel("$Active$");plt.ylabel("$Passive$")
    plt.imshow(param_space,cmap='jet',origin='lower',interpolation='none')
    plt.title("$Core\ Size:\ {0}$".format(core))
    plt.colorbar()
    plt.show()
    return


if __name__=='__main__':
    T=50
    active=6; passive=17 
    cycle=active+passive; num_states=cycle+1
    core_len = 2; core_size = np.square(core_len)
##    emb_core_len = core_len+2; emb_core_size=np.square(emb_core_len)
    
    cmap = colors.ListedColormap(['xkcd:pale grey','xkcd:darkish red','xkcd:almost black'])
    bounds = [0,0.99,active+0.99,cycle+0.99]
    norm = colors.BoundaryNorm(bounds,cmap.N)

    data_path = "./result/"

    config_ids = np.load(data_path+"act_config_ids_states-({0:02d},{1:02d})_core-{2:02d}.npy".format(active,passive,core_size))
    print('Total configs =',num_states**core_size,'\nActive configs = ',len(config_ids))
    print("Persistent config IDs for (active,passive) = ({0:02d},{1:02d}) and core size = {2:02d}\n".format(active,passive,core_size),config_ids)
    
    conf_id = 1496#config_ids[73]
    s,p = run_config_id((active,passive),core_size,conf_id,T=T)#,emb=True,return_str=False)
    com = [ca.Population(p=p,act=active,pas=passive)]
    ca.plot(com,cid=conf_id,cbar='v',txt='True')
    ca.animate(com,s,T=T,interval=500,txt=True)

##    flag='Fixed\ active';anch = 17
##    plot_perst_prob(core_size,flag,anch)
    for i in range(1,18):
        flag='Fixed\ active'
##        flag='Fixed\ passive'
        plot_perst_prob(core_size,flag,i,show=False)
        plt.legend(loc='best',title="${}$".format(flag))
    plt.show()
    plot_param_space(core_size)
