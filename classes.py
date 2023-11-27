# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 10:03:12 2022

@author: Adrien Ange Andre Laurent
"""

# Packages
import itertools as it
import math




# Auxiliary functions
def perm_parity(lsta):
    '''\
    Given a permutation of the digits 0..N in order as a tuple,
    returns its parity (or sign): +1 for even parity; -1 for odd.
    '''
    lst=list(lsta)
    parity = 1
    for i in range(0,len(lst)-1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i,len(lst)), key=lst.__getitem__)
            lst[i],lst[mn] = lst[mn],lst[i]
    return parity


def perm_inv(perm):
    inverse = [0] * len(perm)
    for i, p in enumerate(perm):
        inverse[p] = i
    return tuple(inverse)




# Python class for a single aromatic forest
class aromatic_forest:
    def __init__(self, dic=None):
        if dic is None:
            dic = {}
        self.dic = dic

    def copy(self):
        return aromatic_forest(self.dic.copy())

    def getvertices(self):
        res_root=[0 for i in range(self.n_roots())]
        res_vert=[]
        res_covert=[]
        for v in self.dic.keys():
            if self.dic[v][0]=='r':
                res_root[self.dic[v][1]-1]=v
            if v>0:
                res_vert.append(v)
            else:
                res_covert.append(v)
        return [res_root,res_vert,res_covert]

    def n_roots(self):
        n=0
        for v in self.dic:
            if self.dic[v][0]=='r':
                n+=1
        return n

    def n_vertices(self):
        n=0
        for v in self.dic:
            if v>0:
                n+=1
        return n

    def n_covertices(self):
        n=0
        for v in self.dic:
            if v<0:
                n+=1
        return n

    def order(self):
        return len(self.dic)

    def predecessors(self, v):
        set = []
        for u in self.dic:
            if self.dic[u][1] == v and self.dic[u][0]!='r':
                set.append(u)
        return set

    def roots_cycles(self):
        N=self.order()
        List_cycles=[]
        List_roots=[0 for j in range(0,self.n_roots())]
        list_roots_aux=[]
        list_cycles_aux=[]
        for v in self.dic:
            if not(v in list_roots_aux+list_cycles_aux):
                if self.dic[v][0]=='r':
                    List_roots[self.dic[v][1]-1]=v
                    list_roots_aux.append(v)
                else:
                    i=1
                    bool=1
                    l=[]
                    u=v
                    while self.dic[u][0]!='r' and i<=N and bool:
                        u=self.dic[u][1]
                        if u==v:
                            bool=0
                        i+=1
                        l.append(u)
                    if not(bool):
                        List_cycles.append(l)
                        list_cycles_aux+=l
        return [List_cycles,List_roots]


    def draw_predecessors(self,v,cycle):
        print('[',end='')
        list_pred=self.predecessors(v)
        i=1
        for w in list_pred:
            if not(w in cycle):
                if i:
                    self.draw_predecessors(w,[])
                    i=0
                else:
                    print(',',end='')
                    self.draw_predecessors(w,[])
        print(']',end='')


    def draw(self):
        # Implemented only on the bottom line of the bicomplex!!!
        list_rc=self.roots_cycles()
        for cycle in list_rc[0]:
            print('(',end='')
            i=1
            for v in cycle:
                if i:
                    self.draw_predecessors(v,cycle)
                    i=0
                else:
                    print(',',end='')
                    self.draw_predecessors(v,cycle)
            print(') ',end='')
        for r in list_rc[1]:
            self.draw_predecessors(r,[])
            print(' ',end='')


    def root_wedge(self):
        res=aromatic_form()
        Lroots=self.getvertices()[0]
        for sig in it.permutations(range(0,len(Lroots))):
            aux=self.copy()
            for i in range(0,len(Lroots)):
                aux.dic[Lroots[i]]=['r',1+sig[i]]
            res.add_one_elmt(aux,perm_parity(sig))
        res.ext_mult(1/math.factorial(self.n_roots()))
        return res

    def covertices_wedge(self):
        res=aromatic_form()
        for sig in it.permutations(range(0,self.n_covertices())):
            aux=self.copy()
            for v in self.dic:
                if v<0:
                    aux.dic[v]=self.dic[-sig[-v-1]-1]
                if aux.dic[v][1]<0:
                    aux.dic[v]=[aux.dic[v][0],-perm_inv(sig)[-aux.dic[v][1]-1]-1]
            res.add_one_elmt(aux,perm_parity(sig))
            res.simplify()
        res.ext_mult(1/math.factorial(self.n_covertices()))
        return res

    def wedge(self):
        g=self.covertices_wedge().root_wedge()
        g.simplify()
        return g

    def is_equal(self, aro2):
        dic1 = self.dic
        dic2 = aro2.dic
        N=self.n_vertices()
        res = False
        for sig in list(it.permutations(range(0, N))):
            dicaux = {}
            for v in dic1:
                if v>0:
                    dicaux[v]=dic1[sig[v-1]+1]
                else:
                    dicaux[v]=dic1[v]
                if dicaux[v][1]>0 and dicaux[v][0]!='r':
                    dicaux[v]=[dicaux[v][0],perm_inv(sig)[dicaux[v][1]-1]+1]
            if dicaux == dic2:
                res = True
                break
        return res

    def d_H(self):
        Lroots=self.getvertices()[0]
        r=Lroots[-1]
        res=aromatic_form()
        for v in self.dic:
            aux=self.copy()
            aux.dic[r]=['v',v]
            res.add_one_elmt(aux,1)
        res.simplify()
        return res

    def h_H(self):
        res=aromatic_form()
        N=self.order()
        n=self.n_roots()
        for p in range(1,N+1):
            for v in self.dic:
                v_pred=self.predecessors(v)
                pi_v=len(v_pred)
                if pi_v>=p:
                    list_nodes=[w for w in self.dic]
                    list_other_nodes=[w for w in self.dic]
                    list_other_nodes.remove(v)
                    for r0 in v_pred:
                        for sublist in list(map(list, it.combinations([x for x in v_pred if x!=r0], pi_v-p))):
                            comp_sublist=[item for item in v_pred if item not in sublist and item!=r0]
                            if pi_v-p==1:
                                graftings=[[item] for item in list_other_nodes]
                            else:
                                graftings=[item for item in it.product(list_other_nodes,repeat=pi_v-p)]
                            for f in graftings:
                                if p==2:
                                    graftings2=[[item] for item in list_nodes]
                                else:
                                    graftings2=[item for item in it.product(list_nodes,repeat=p-1)]
                                for f2 in graftings2:
                                    gaux=self.copy()
                                    for j in range(0,len(f)):
                                        gaux.dic[sublist[j]]=['v',f[j]]
                                    for j in range(0,len(f2)):
                                        gaux.dic[comp_sublist[j]]=['v',f2[j]]
                                    gaux.dic[r0]=['r',n+1]
                                    res.add_one_elmt(gaux,(-1)**(pi_v-p)/(p+n)*(n+1))
        res.ext_mult(1/N)
        res.simplify()
        return res.wedge()

    def d_V(self):
        p=self.n_covertices()
        res=aromatic_form()
        for v in self.dic:
            if v>0:
                aux=aromatic_forest()
                for w in self.dic:
                    if w!=v:
                        a=self.dic[w]
                        boo=(w>v)
                        if a[1]>v and a[0]!='r':
                            aux.dic[w-boo]=['v',a[1]-1]
                        else:
                            if a[1]==v and a[0]!='r':
                                aux.dic[w-boo]=['v',-p-1]
                            else:
                                aux.dic[w-boo]=[a[0],a[1]]
                b=self.dic[v]
                if b[1]>v and b[0]!='r':
                    aux.dic[-p-1]=['v',b[1]-1]
                else:
                    if b[1]==v and b[0]!='r':
                        aux.dic[-p-1]=['v',-p-1]
                    else:
                        aux.dic[-p-1]=[b[0],b[1]]
                res.add_one_elmt(aux,1)
        res.simplify()
        return res.covertices_wedge()

    def h_V(self):
        p=self.n_covertices()
        m=self.n_vertices()
        aux=aromatic_forest()
        for v in self.dic:
            if v!=-p:
                a=self.dic[v]
                if a[1]==-p:
                    aux.dic[v]=[a[0],m+1]
                else:
                    aux.dic[v]=[a[0],a[1]]
        b=self.dic[-p]
        if b[1]==-p:
            aux.dic[m+1]=[b[0],m+1]
        else:
            aux.dic[m+1]=[b[0],b[1]]
        res=aromatic_form([[aux,1]])
        res.ext_mult(p/self.order())
        return res

    def Euler_op(self):
        res=aromatic_form()
        for v in self.dic:
            v_pred=self.predecessors(v)
            list_other_nodes=[w for w in self.dic]
            list_other_nodes.remove(v)
            if len(v_pred)==1:
                graftings=[[item] for item in list_other_nodes]
            else:
                graftings=[item for item in it.product(list_other_nodes,repeat=len(v_pred))]
            for f in graftings:
                gaux=self.copy()
                for j in range(0,len(f)):
                    gaux.dic[v_pred[j]]=['v',f[j]]
                gaux.dic[v]=[gaux.dic[v][0],gaux.dic[v][1]]
                res.add_one_elmt(gaux,(-1)**(len(v_pred)))
        res.simplify()
        return res

    def Euler_op_O(self):
        res=aromatic_form()
        N=self.order()
        for v in self.dic:
            list_pred=self.predecessors(v)
            v_pred=[]
            for w in list_pred:
                if w<v:
                    v_pred.append(w)
                else:
                    if w==v:
                        v_pred.append(-1)
                    else:
                        v_pred.append(w-1)
            list_other_nodes=range(1,N)
            if len(v_pred)==1:
                graftings=[[item] for item in list_other_nodes]
            else:
                graftings=[item for item in it.product(list_other_nodes,repeat=len(v_pred))]
            for f in graftings:
                gaux=aromatic_forest()
                for j in range(0,len(f)):
                    gaux.dic[v_pred[j]]=['v',f[j]]
                for w in list_other_nodes:
                    if not(w in v_pred):
                        x=self.dic[w+(w>=v)][1]
                        gaux.dic[w]=['v',x-(x>v)]
                if not(-1 in v_pred):
                    x=self.dic[v][1]
                    gaux.dic[-1]=['v',x-(x>v)]
                res.add_one_elmt(gaux,(-1)**(len(v_pred)))
        res.simplify()
        return res

    def decomp_higher_Euler_op(self):
        res=aromatic_form()
        N=self.order()
        for p in range(0,N+1):
            for v in self.dic:
                v_pred=self.predecessors(v)
                pi_v=len(v_pred)
                if pi_v>=p:
                    list_nodes=[w for w in self.dic]
                    list_other_nodes=[w for w in self.dic]
                    list_other_nodes.remove(v)
                    for sublist in list(map(list, it.combinations(v_pred, pi_v-p))):
                        comp_sublist=[item for item in v_pred if item not in sublist]
                        if pi_v-p==1:
                            graftings=[[item] for item in list_other_nodes]
                        else:
                            graftings=[item for item in it.product(list_other_nodes,repeat=pi_v-p)]
                        for f in graftings:
                            if p==1:
                                graftings2=[[item] for item in list_nodes]
                            else:
                                graftings2=[item for item in it.product(list_nodes,repeat=p)]
                            for f2 in graftings2:
                                gaux=self.copy()
                                for j in range(0,len(f)):
                                    gaux.dic[sublist[j]]=['v',f[j]]
                                for j in range(0,len(f2)):
                                    gaux.dic[comp_sublist[j]]=['v',f2[j]]
                                res.add_one_elmt(gaux,(-1)**(pi_v-p))
        res.ext_mult(1/N)
        res.simplify()
        return res

    def I(self):
        res=aromatic_form()
        p=self.n_covertices()
        v_pred=self.predecessors(-p)
        list_other_nodes=[w for w in self.dic]
        list_other_nodes.remove(-p)
        if len(v_pred)==1:
            graftings=[[item] for item in list_other_nodes]
        else:
            graftings=[item for item in it.product(list_other_nodes,repeat=len(v_pred))]
        for f in graftings:
            gaux=self.copy()
            for j in range(0,len(f)):
                gaux.dic[v_pred[j]]=['v',f[j]]
            gaux.dic[-p]=['v',gaux.dic[-p][1]]
            res.add_one_elmt(gaux,(-1)**(len(v_pred)))
        res.simplify()
        return res.covertices_wedge()

    def h_I(self):
        res=aromatic_form()
        p=self.n_covertices()
        v_pred=self.predecessors(-p)
        pi_v=len(v_pred)
        for q in range(1,pi_v+1):
            list_nodes=[w for w in self.dic]
            list_other_nodes=[w for w in self.dic]
            list_other_nodes.remove(-p)
            for r0 in v_pred:
                for sublist in list(map(list, it.combinations([x for x in v_pred if x!=r0], pi_v-q))):
                    comp_sublist=[item for item in v_pred if item not in sublist and item!=r0]
                    if pi_v-q==1:
                        graftings=[[item] for item in list_other_nodes]
                    else:
                        graftings=[item for item in it.product(list_other_nodes,repeat=pi_v-q)]
                    for f in graftings:
                        if q==2:
                            graftings2=[[item] for item in list_nodes]
                        else:
                            graftings2=[item for item in it.product(list_nodes,repeat=q-1)]
                        for f2 in graftings2:
                            gaux=self.copy()
                            for j in range(0,len(f)):
                                gaux.dic[sublist[j]]=['v',f[j]]
                            for j in range(0,len(f2)):
                                gaux.dic[comp_sublist[j]]=['v',f2[j]]
                            gaux.dic[r0]=['r',1]
                            res.add_one_elmt(gaux,(-1)**(pi_v-q)/q)
        res.simplify()
        return res.covertices_wedge()

    def delta_V(self):
        return self.d_V().I()

    def h_delta_V(self):
        return self.h_V().I()

    def is_1_loop(self):
        aro_boo=False
        for v in self.dic:
            if self.dic[v][0]=='v' and self.dic[v][1]==v:
                aro_boo=True
                break
        return aro_boo
    
    def is_linear(self):
        aro_boo=True
        for v in self.dic:
            if len(self.predecessors(v))>1:
                aro_boo=False
                break
        return aro_boo

    def height(self):
        gdic=self.dic
        res=[]
        N=len(gdic)
        for v in range(1,N+1):
            h=0
            w=v
            l=[v]
            while gdic[w][0]!='r':
                h+=1
                w=gdic[w][1]
                if w in l:
                    h=-1
                    break
                else:
                    l.append(w)
            res.append(h)
        return res

    def h_H_div_free(self):
        res=aromatic_form()
        N=self.order()
        n=self.n_roots()
        for p in range(1,N+1):
            for v in self.dic:
                v_pred=self.predecessors(v)
                pi_v=len(v_pred)
                if pi_v>=p:
                    list_nodes=[w for w in self.dic]
                    list_other_nodes=[w for w in self.dic]
                    list_other_nodes.remove(v)
                    for r0 in v_pred:
                        for sublist in list(map(list, it.combinations([x for x in v_pred if x!=r0], pi_v-p))):
                            comp_sublist=[item for item in v_pred if item not in sublist and item!=r0]
                            if pi_v-p==1:
                                graftings=[[item] for item in list_other_nodes]
                            else:
                                graftings=[item for item in it.product(list_other_nodes,repeat=pi_v-p)]
                            for f in graftings:
                                if p==2:
                                    graftings2=[[item] for item in list_nodes]
                                else:
                                    graftings2=[item for item in it.product(list_nodes,repeat=p-1)]
                                for f2 in graftings2:
                                    gaux=self.copy()
                                    for j in range(0,len(f)):
                                        gaux.dic[sublist[j]]=['v',f[j]]
                                    for j in range(0,len(f2)):
                                        gaux.dic[comp_sublist[j]]=['v',f2[j]]
                                    gaux.dic[r0]=['r',n+1]
                                    if self.dic[v][0]=='r':
                                        res.add_one_elmt(gaux,(-1)**(pi_v-p)/(p+n-1)*(n+1))
                                    else:
                                        res.add_one_elmt(gaux,(-1)**(pi_v-p)/(p+n)*(n+1))
        res.ext_mult(1/N)
        res.simplify()
        res.simplify_1_loop()
        return res.wedge()

    def correction_horizontal_homotopy_div_free(self):
        '''
        Returns the additional term L_r/N in the horizontal homotopy in the divergence-free case with n=1 root.
        '''
        res=aromatic_form()
        N=self.order()
        n=self.n_roots()
        if n==1:
            pmax=0
        else:
            pmax=-1
        for p in range(0,pmax+1):
            for v in self.getvertices()[0]:
                v_pred=self.predecessors(v)
                pi_v=len(v_pred)
                if pi_v>=p:
                    list_nodes=[w for w in self.dic]
                    list_other_nodes=[w for w in self.dic]
                    list_other_nodes.remove(v)
                    for sublist in list(map(list, it.combinations(v_pred, pi_v-p))):
                        comp_sublist=[item for item in v_pred if item not in sublist]
                        if pi_v-p==1:
                            graftings=[[item] for item in list_other_nodes]
                        else:
                            graftings=[item for item in it.product(list_other_nodes,repeat=pi_v-p)]
                        for f in graftings:
                            if p==1:
                                graftings2=[[item] for item in list_nodes]
                            else:
                                graftings2=[item for item in it.product(list_nodes,repeat=p)]
                            for f2 in graftings2:
                                gaux=self.copy()
                                for j in range(0,len(f)):
                                    gaux.dic[sublist[j]]=['v',f[j]]
                                for j in range(0,len(f2)):
                                    gaux.dic[comp_sublist[j]]=['v',f2[j]]
                                res.add_one_elmt(gaux,(-1)**(pi_v-p))
        res.ext_mult(1/N)
        res.simplify()
        res.simplify_1_loop()
        return res

    def Euler_op_r_div_free(self):
        res=aromatic_form()
        v=self.getvertices()[0][0]
        v_pred=self.predecessors(v)
        list_other_nodes=[w for w in self.dic]
        list_other_nodes.remove(v)
        if len(v_pred)==1:
            graftings=[[item] for item in list_other_nodes]
        else:
            graftings=[item for item in it.product(list_other_nodes,repeat=len(v_pred))]
        for f in graftings:
            gaux=self.copy()
            for j in range(0,len(f)):
                gaux.dic[v_pred[j]]=['v',f[j]]
            gaux.dic[v]=[gaux.dic[v][0],gaux.dic[v][1]]
            res.add_one_elmt(gaux,(-1)**(len(v_pred)))
        res.simplify()
        res.simplify_1_loop()
        return res

    def Euler_op_O_r_div_free(self):
        res=aromatic_form()
        N=self.order()
        v=self.getvertices()[0][0]
        list_pred=self.predecessors(v)
        v_pred=[]
        for w in list_pred:
            if w<v:
                v_pred.append(w)
            else:
                if w==v:
                    v_pred.append(-1)
                else:
                    v_pred.append(w-1)
        list_other_nodes=range(1,N)
        if len(v_pred)==1:
            graftings=[[item] for item in list_other_nodes]
        else:
            graftings=[item for item in it.product(list_other_nodes,repeat=len(v_pred))]
        for f in graftings:
            gaux=aromatic_forest()
            for j in range(0,len(f)):
                gaux.dic[v_pred[j]]=['v',f[j]]
            for w in list_other_nodes:
                if not(w in v_pred):
                    x=self.dic[w+(w>=v)][1]
                    gaux.dic[w]=['v',x-(x>v)]
            if not(-1 in v_pred):
                x=self.dic[v][1]
                gaux.dic[-1]=['v',x-(x>v)]
            res.add_one_elmt(gaux,(-1)**(len(v_pred)))
        res.simplify()
        res.simplify_1_loop()
        return res

    def alternative_h_H(self):
        aux=self.Euler_op()
        aux.ext_mult(-1/self.order())
        aux.add_one_elmt(self,1)
        res=aromatic_form()
        while aux.is_1_loop():
            for x in aux.alist:
                if x[0].is_1_loop():
                    for v in x[0].dic:
                        if x[0].dic[v][0]=='v' and x[0].dic[v][1]==v:
                            y=aromatic_forest()
                            for w in x[0].dic:
                                if w==v:
                                    y.dic[v]=['r',1]
                                else:
                                    y.dic[w]=x[0].dic[w]
                            res.add_one_elmt(y,x[1])
                            corr=y.d_H()
                            corr.ext_mult(-x[1])
                            aux.add(corr)
                            break
                    break
        res.simplify()
        return res

    def i_Lie(self,vf):
        p=self.n_covertices()
        m=self.n_vertices()
        res=aromatic_form()
        if p==0:
            return res
        else:
            for vfformtot in vf.alist:
                vfform=vfformtot[0]
                aux=aromatic_forest()
                liste_graft=[]
                for v in vfform.dic:
                    if v!=1:
                        aux.dic[m+v]=['v',vfform.dic[v][1]+m]
                    else:
                        b=self.dic[-p]
                        if b[1]==-p:
                            aux.dic[m+1]=[b[0],0]
                            liste_graft.append(m+1)
                        else:
                            aux.dic[m+1]=[b[0],b[1]]
                for v in self.dic:
                    if v!=-p:
                        a=self.dic[v]
                        if a[1]==-p:
                            aux.dic[v]=[a[0],0]
                            liste_graft.append(v)
                        else:
                            aux.dic[v]=[a[0],a[1]]
                liste_target=[m+i for i in range(1,1+vfform.order())]
                if len(liste_graft)==1:
                    graftings=[[item] for item in liste_target]
                else:
                    graftings=[item for item in it.product(liste_target,repeat=len(liste_graft))]
                for f in graftings:
                    gaux=aux.copy()
                    for j in range(0,len(f)):
                        gaux.dic[liste_graft[j]]=['v',f[j]]
                    res.add_one_elmt(gaux,vfformtot[1]*p)
            res.simplify()
            return res

    def Lie(self,vf):
        aux=self.d_V().i_Lie(vf)
        aux.add(self.i_Lie(vf).d_V())
        aux.simplify()
        return aux

    def h_H_Lie(self,vf):
        return self.d_V().h_H().i_Lie(vf)

    def h_I_Lie(self,vf):
        return self.d_V().h_I().i_Lie(vf)



















# Python class for aromatic forms
class aromatic_form:
    def __init__(self, alist=None):
        if alist is None:
            alist = []
        self.alist = alist

    def add_one_elmt(self, aromatic_forest, cst):
        self.alist.append([aromatic_forest, cst])

    def freelist(self):
        res = []
        for aro in self.alist:
            res.append([aro[0].dic, aro[1]])
        return res

    def aro_print(self):
        aux = self.freelist()
        print(*aux, sep="\n")

    def draw(self):
        for aro in self.alist:
            print(round(aro[1],5),end=' ')
            aro[0].draw()
            print()

    def ext_mult(self, c):
        if c==0:
            self.alist=[]
        else:
            for g in self.alist:
                g[1] = c*g[1]

    def add(self, alistaux):
        self.alist+=alistaux.alist
        self.simplify()

    def add_no_simplify(self, alistaux):
        self.alist+=alistaux.alist

    def simplify(self):
        ll=len(self.alist)
        i = 0
        while i<ll-1:
            j=i+1
            while j<ll:
                if self.alist[i][0].is_equal(self.alist[j][0]):
                    self.alist[i][1]+=self.alist[j][1]
                    del self.alist[j]
                    ll-=1
                else:
                    j+=1
            if abs(self.alist[i][1])<10**(-10):
                del self.alist[i]
                ll-=1
            else:
                i+=1
        if ll==1 and abs(self.alist[0][1])<10**(-10):
            del self.alist[0]

    def root_wedge(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].root_wedge()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def covertices_wedge(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].covertices_wedge()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def wedge(self):
        g=self.covertices_wedge().root_wedge()
        return g

    def d_H(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].d_H()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def d_V(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].d_V()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def h_V(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].h_V()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def Euler_op(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].Euler_op()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def Euler_op_O(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].Euler_op_O()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def h_H(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].h_H()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def decomp_higher_Euler_op(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].decomp_higher_Euler_op()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def I(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].I()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def h_I(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].h_I()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def delta_V(self):
        return self.d_V().I()

    def h_delta_V(self):
        return self.h_V().I()
    
    def is_1_loop(self):
        boo=False
        for x in self.alist:
            if x[0].is_1_loop():
                boo=True
                break
        return boo

    def simplify_1_loop(self):
        ll=len(self.alist)
        i = 0
        while i<ll:
            if self.alist[i][0].is_1_loop():
                del self.alist[i]
                ll-=1
            else:
                i+=1
    
    def simplify_linear(self):
        ll=len(self.alist)
        i = 0
        while i<ll:
            if self.alist[i][0].is_linear():
                i+=1
            else:
                del self.alist[i]
                ll-=1

    def d_H_div_free(self):
        res=self.d_H()
        res.simplify_1_loop()
        return res

    def d_V_div_free(self):
        res=self.d_V()
        res.simplify_1_loop()
        return res

    def h_V_div_free(self):
        res=self.h_V()
        res.simplify_1_loop()
        return res

    def Euler_op_div_free(self):
        res=self.Euler_op()
        res.simplify_1_loop()
        return res

    def I_div_free(self):
        res=self.I()
        res.simplify_1_loop()
        return res

    def h_I_div_free(self):
        res=self.h_I()
        res.simplify_1_loop()
        return res

    def delta_V_div_free(self):
        return self.d_V_div_free().I_div_free()

    def h_delta_V_div_free(self):
        return self.h_V_div_free().I_div_free()

    def h_H_div_free(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].h_H_div_free()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def correction_horizontal_homotopy_div_free(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].correction_horizontal_homotopy_div_free()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def Euler_op_r_div_free(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].Euler_op_r_div_free()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def Euler_op_O_r_div_free(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].Euler_op_O_r_div_free()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def alternative_h_H(self):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].alternative_h_H()
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def i_Lie(self,vf):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].i_Lie(vf)
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def Lie(self,vf):
        res=aromatic_form()
        for g in self.alist:
            aux=g[0].Lie(vf)
            aux.ext_mult(g[1])
            res.add_no_simplify(aux)
        res.simplify()
        return res

    def h_H_Lie(self,vf):
        return self.d_V().h_H().i_Lie(vf)

    def h_I_Lie(self,vf):
        return self.d_V().h_I().i_Lie(vf)









