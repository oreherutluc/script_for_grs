# -*- coding:utf-8 -*-
from math import cos,sin,sqrt,pi
from copy import deepcopy

cos_deg=lambda x: cos(x*pi/180.0)
sin_deg=lambda x: sin(x*pi/180.0)

class GULP():
    def __init__(self, filename=None):
        #открывается файл
        f=open(filename,'r')
        x=f.readlines()
        f.close()
        gulp_list=[]
        i=0
        while True:
            elem=x[i]
            if elem.strip()=='cell':
                i+=1
                # cell_str -- это объект класса cell, содержащий в своих ПОЛЯХ информацию о ячейке из файла
                cell_str = cell(x[i])
            # gulp_list -- это список объектов класса gulp_atom, содержащих в своих ПОЛЯХ информацию о каждом атоме файла
            elif elem[0:10].strip()=='fractional':
                i+=1
                while x[i].strip()!='space':
                    gulp_atom_str=gulp_atom(x[i])
                    gulp_list.append(gulp_atom_str)
                    i+=1
            i+=1
            if i==len(x):
                break
        self.cellmethod = [cell_str.a, cell_str.b, cell_str.c, cell_str.alfa, cell_str.beta, cell_str.gamma]
        #в поле cell записывается объект класса cell
        self.cell = cell_str
        #в поле atoms объекта класса GULP записывается список объектов класса gulp_atom
        self.atoms = gulp_list

    def gulp2xyz(self):
        A11=self.cell.a
        A21=self.cell.b*cos_deg(self.cell.gamma)
        A22=self.cell.b*sin_deg(self.cell.gamma)
        A31=self.cell.c*cos_deg(self.cell.beta)
        A32=self.cell.c*((cos_deg(self.cell.alfa)-cos_deg(self.cell.beta)*cos_deg(self.cell.gamma))/sin_deg(self.cell.gamma))
        fir=self.cell.c*self.cell.a*self.cell.b
        sec=1.0-(cos_deg(self.cell.alfa))**2-(cos_deg(self.cell.beta))**2-(cos_deg(self.cell.gamma))**2
        forth=2.0*cos_deg(self.cell.alfa)*cos_deg(self.cell.beta)*cos_deg(self.cell.gamma)
        omega=fir*sqrt(sec+forth)
        A33=omega/(self.cell.a*self.cell.b*sin_deg(self.cell.gamma))
        xyzs=[]
        #elem -- это объект класса gulp_atom, описывающий один атом
        for elem in self.atoms:
            # применяем метод make_xyz, который на вход получает объект класса gulp_atom и объект класса cell и выводит
            # декартовы координаты одного атома в виде отформатированной строки
            # xyz - перменная, в которую записана строка с декартовыми координатами одного атома
            try:
                Dx=A11*elem.x
                Dy=A21*elem.x+A22*elem.y
                Dz=A31*elem.x+A32*elem.y+A33*elem.z
            except AttributeError:
                print
            xyzs+=['{0:2s} {1:11.9f} {2:11.9f} {3:11.9f}'.format(elem.elem, Dx, Dy, Dz)]
            #xyz=elem.make_xyz(self.cell)
            #xyzs.append(xyz)
        return xyzs

#без __add__ нельзя применить + к двум объектам
    def __add__(self, other):
        #новый экземпляр копии self, вторая ссылка на него
        total=deepcopy(self)
        total.atoms+=other.atoms
        return total

    def droptofile(self, filename=None):
        if filename!=None:
            out=str(self.cell)+'\nfractional\n'+'\n'.join(map(str,self.atoms))
            dfile=open(filename, 'w')
            dfile.write(out)
            dfile.close()
#map(str,self.atoms)

class gulp_atom():
    @staticmethod
    def stroka(s):
        res=s.split()
        for i in xrange(len(res)):
            try:
                res[i]=float(res[i])
            except ValueError:
                pass

        return res
    def __init__(self, line=None):
        if line==None:
            self.elem=''
            self.elem_class=''
            self.x=0.0
            self.y=0.0
            self.z=0.0
            self.charge=0.0
            self.one=0.0
            self.zero=0.0
        else:
            spstroka=line.split()
            if spstroka[0].strip() != 'space' and len(spstroka)==8 :
                self.elem, self.elem_class, self.x, self.y, self.z, self.charge, self.one, self.zero=self.stroka(line)
                # try: self.x, self.y, self.z, self.charge, self.one, self.zero=map(float, spstroka[2:])
                # except ValueError:
                #     print('В строке неверные данные {0:s}'.format(line))
                #     raise SystemExit(1)
    def __str__(self):
        return '{0:2s} {1:6s} {2:11.9f} {3:11.9f} {4:11.9f} {5:11.9f} {6:7.5f} {7:7.5f}'.format(self.elem, self.elem_class, self.x, self.y, self.z, self.charge, self.one, self.zero)

class cell():
    def __init__(self, line=None):
        self.line = None
        if line==None:
            self.a=0.0
            self.b=0.0
            self.c=0.0
            self.alfa=0.0
            self.beta=0.0
            self.gamma=0.0
        else:
            q=line.split()
            if len(q)==6:
                try:
                    self.a, self.b, self.c, self.alfa, self.beta, self.gamma=map(float,q)
                except ValueError:
                    print 'в строке не те данные {0:s}'.format(line)
                    raise SystemExit(1)
            else:
                print 'в строке не хватает данных: {0:s}'.format(line)
                raise SystemExit(1)
        self.space=1
    def __str__(self):
        return '{0:4s} \n {1:8.4f} {2:8.4f} {3:8.4f} {4:8.3f} {5:8.3f} {6:8.3f}'.format('cell', self.a, self.b, self.c, self.alfa, self.beta, self.gamma)


testabc=cell(line='9.5108 9.5108 12.991 90 90 120')
testatom=gulp_atom(line='Si        core      0.115329999  0.730440000  0.381549999  2.3999999  1.0000  0.0000')
g=GULP(filename='K2Si2O5_geometry_simul_fit_s_1_1.grs')
m=g.gulp2xyz()
#print m

t=cell(line='9.5108 9.5108 12.991 90 90 120')
print(testabc.__str__())
x=g+GULP(filename='Al2O3_young_relaxed_fit_s_2_2.grs')
print x

g.droptofile(filename='tmp.grs')
#z.write_to_xyz('file.xyz')
