import numpy as np
import ghca_main as ca

def base_conv(number,base=2,padding=0):
    digits = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    num = abs(number)
    res = []
    while num:
        res.append(digits[int(num % base)])
        num //= base
    if padding:
        res.append('0' * padding)
    if number < 0:
        res.append('-')
    return ''.join(reversed(res or '0'))

def str_to_state(str_config,base):
    pop = []
    for s in str_config:
        pop.append(int(s,base))
    pop = np.array(pop,dtype=np.int8)
    pop = np.reshape(pop,(int(np.sqrt(len(str_config))),int(np.sqrt(len(str_config)))))
    return pop

def state_to_str(pop,base):
    str_config = ''
    pop = pop.flatten()
    for i in pop:
        str_config+=base_conv(i,base)
    return str_config

def emb(ic,core_size,active,passive,T):
    n_states = active+passive+1
    grid_size = np.sqrt(core_size) + 2
    s_str = base_conv(ic,n_states);s_str=s_str.rjust(core_size,'0')
    q = str_to_state(s_str,n_states)
    p = np.zeros((int(grid_size),int(grid_size)),dtype=np.int8)
    if int(np.sqrt(core_size))%2==0:
        p[int(grid_size/2)-int(np.sqrt(core_size)/2):int(grid_size/2)+int(np.sqrt(core_size)/2),
        int(grid_size/2)-int(np.sqrt(core_size)/2):int(grid_size/2)+int(np.sqrt(core_size)/2)]=q
    else:
        p[int(grid_size/2)-int(np.sqrt(core_size)/2):int(grid_size/2)+int(np.sqrt(core_size)/2)+1,
        int(grid_size/2)-int(np.sqrt(core_size)/2):int(grid_size/2)+int(np.sqrt(core_size)/2)+1]=q

    com = [ca.Population(p=np.copy(p),act=active,pas=passive,periodic=True)]
    states = np.zeros((T,int(grid_size),int(grid_size)),dtype=np.int8)
    for t in range(T):
        states[t] = np.copy(com[0].p)
        com = ca.run(com)
    ca.animate(com,states,T=T,txt=True,interval=500)
    return

if __name__=='__main__':
    pass
