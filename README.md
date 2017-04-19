# proteinfolding

主程序IDPmd13_2.f90

进行动力学模拟的输入文件：
1.构象坐标文件（文件名N-PDBid.dat），包含每个氨基酸残基C-alpha原子的三维坐标； 2.相互作用文件（native contact map，文件名appNCS_PDBid.dat），包含每队残基的相互作用倾向和天然结构中的距离；
3.参数文件（input.dat）

参数文件使用说明见word文档。

以上输入文件用于蛋白质复合物的解离过程，每个坐标均包含两个蛋白质单体。
附件中的五个体系大小不同，可以用于测试加速效果。

如果需要模拟其它体系，需要先用PDB_Contact_Map程序对原始PDB文件进行处理。

祝使用愉快

配置文件说明:
===============================================================================
1STF      10-Apr-2017    
.1STF.test    
95  307  20 1    
N-1STF.dat    
N-1STF.dat    
appNCS-1STF.dat    
-----------------parameter--------------------    
1334887077757897                                 ; random numer, 16 characters    
      1.000       1.000       1.000      1.000   ; epsil, epsil1, epsil2, enscale    
    100.000      20.000       1.000      0.500   ; ck_r, ck_tht, ck_phi1, ck_phi3    
      4.000       1.000       1.200              ; sigma_ij, amass, gamma    
      1.000       1.000                          ; avsn0, boltz    
      1.000       0.100                          ; gr0, gdr    
      0.005       1.000       0                  ; s_dt, s_gm, iFixKr    
      6           2           3                  ; k_sol, n_sol, m_sol    
      0.200       0.100                          ; epsilon_p, epsilon_pp    
      1                                          ; temp    
      0.0         0.0         0.0                ; ga1_f1, ga2_f1, gQ0_f1    
      0.0         0.0         0.0                ; ga1_f2, ga2_f2, gQ0_f2    
      0.0         0.0         0.0                ; ga1_b, ga2_b, gQ0_b    
      0.0         0.0      0.0       0.02        ; ga1_w, ga2_w, gQ0_w, alpha_Qw    
      1                                          ; Is_Solvation    
      1           0           200                ; nConform, nCon0, nRunConf    
      500 0.01  0 1000  0 1000                   ; gQbnativ, gQbdenatural, gQf1nativ, gQf1denatural, gQf2nativ, gQf2denatural    
     500000     50000000                         ; nsnap, nstep    
        400      100000    500000                ; nConformOutput, nOutput0, ndOutput    
-2 1000  -2 1000  -2 1000                        ; outQf1_i, outQf1_f, outQf2_i, outQf2_f, outQb_i, outQb_f    
         20     2000000                          ; nbinsnap, nbinsnap0    
        100         100                          ; nbin_f, nbin_b    
      3       1        0.000                     ; dbin_f, dbin_b, vbin0    
          0         100         7.0   -350.0     ; IsEbin, nEbin, dEbin, vEbin0    
          1          100         0.6   -57       ; IsEbbin, nEbbin, dEbbin, vEbbin0    
          0          100        1.000     0.000  ; IsRbin, nRbin, dRbin, vRbin0    
          1          100        0.5   -25  0.001 ; IsWbin, nWbin, dWbin, vWbin0, cri_Qb    
      7500.0      10.000                         ; PBC: pL, dl    
      1.200       1.200         0.96             ; Alpha1, Alpha2, Beta    
      0.000       5.000                          ; Delta, CritR_non    
参数说明：    
    
1STF      10-Apr-2017        
模拟的体系名与模拟日期。不影响程序的运行。    
.1STF.test    
统一的后缀名。方便整理用，本次模拟所有的输出文件均具有此后缀。    
95  307  20 1    
1.    前两个参数分别是第一条链的长度和总长度，若模拟单条链则两个值相等。    
2.    第三个参数是相互作用文件中用于区分两条链的间隔值。因为相互作用文件中只有残基编号而没有链编号，所以对残基进行了重新编号，并且让两条链的编号中间存在一定的间隔，以便在程序中能够通过读入残基编号来判断其属于哪一条链。以该体系为例，编号1-95的残基为A链，编号116以后为B链。建议不修改此参数。    
3.    模拟模式。该版本的程序有0和1两个赋值，模式0为单链模拟，1为双链模拟。    
N-1STF.dat    
天然构象的坐标文件。    
N-1STF.dat    
初始构象的坐标文件。本输入文件用于模拟蛋白质复合物解离的动力学过程，所以初始构象与天然构象相同。模拟开始时的初始构象根据需要而定，不一定与天然构象一致。    
appNCS-1STF.dat    
相互作用文件。    

1334887077757897                                  ; random numer, 16 characters    
随机数种子。    
      1.000       1.000       1.000      1.000    ; epsil, epsil1, epsil2, enscale    
控制力场中模拟强度的参数epsilon。    
    100.000      20.000       1.000      0.500    ; ck_r, ck_tht, ck_phi1, ck_phi3    
力场中的力常数。分别为键能项（r）、键角项（theta）和二面角项（phi）。    
      4.000       1.000       1.200               ; sigma_ij, amass, gamma    
用于朗之万动力学模拟的参数。该模型的动力学过程需要考虑溶液中的布朗运动，因此采用朗之万动力学。    
1.000       1.000                           ; avsn0, boltz    
阿伏伽德罗常数与玻尔兹曼常数。该模型中定义为1。    
参考文献：http://doi.org/10.1016/S0022-2836(02)01434-1    
      1.000       0.100                           ; gr0, gdr    
相互接触项。用于计算残基间接触程度随距离的变化。    
      0.005       1.000       0                   ; s_dt, s_gm, iFixKr    
      6           2           3                   ; k_sol, n_sol, m_sol    
      0.200       0.100                           ; epsilon_p, epsilon_pp    
以上为溶剂化模型参数。    
1    ; temp    
温度。该模型定义为1。参考文献同上。    
    
      0.0         0.0         0.0                 ; ga1_f1, ga2_f1, gQ0_f1    
      0.0         0.0         0.0                 ; ga1_f2, ga2_f2, gQ0_f2    
      0.0         0.0         0.0                 ; ga1_b, ga2_b, gQ0_b    
      0.0         0.0      0.0       0.02         ; ga1_w, ga2_w, gQ0_w, alpha_Qw    
偏置场（bias potential）参数。为了能够充分采样而采用的额外偏置场，对所有偏置场条件模拟完毕之后再将结果整合。用于热力学模拟，因为不能够反映时间所以动力学模拟不采用偏置场。    
      1                                           ; Is_Solvation    
溶剂化模式。0为非溶剂化模型，1为溶剂化模型。    
    
      1           0           200                   ; nConform, nCon0, nRunConf    
本次模拟的构象数、模拟开始前需要跳过的构象数、每个构象的模拟次数。    
      500  0.01  0 1000  0 1000                      ; gQbnativ, gQbdenatural, gQf1nativ, gQf1denatural, gQf2nativ, gQf2denatural    
根据反应坐标Q定义的模拟范围。Q为接触数，两个残基间距小于4.5即认为其存在接触，单链内的接触总数为Qf，两条链之间的接触数为Qb。    
模拟的起止由该项参数决定，体系的Q值在天然态（native）和变性态（denatural）之间时模拟正常进行，在此范围之外则认为已经达到目标状态而终止模拟。根据实际需要设定native和denatural的值。    
     500000     50000000                           ; nsnap, nstep    
Output文件的输出间隔，总模拟步数。    
        400      100000    500000                 ; nConformOutput, nOutput0, ndOutput    
需要输出的构象数目，开始输出构象的步数，构象输出间隔。    
-2 1000  -2 1000  -2 1000                               ;outQf1_i, outQf1_f, outQf2_i, outQf2_f, outQb_i, outQb_f    
构象输出范围。定义类似模拟范围。    
         20     2000000                           ; nbinsnap, nbinsnap0    
采样时间间隔与开始进行采样的步数。    
计算力与速度、采样与信息输出是三个独立的过程。力与速度在每一步都会进行计算，但是只有在到达采样间隔时才会把当前的计算结果输出到相应的文件中。采样的频率直接影响通过输出文件得到的采样分布精度，进而影响模拟结果。信息输出（Output文件）仅仅用于监测模拟是否正常进行，并不用于后续分析。    
        100         100                           ; nbin_f, nbin_b    
每次采样的区间划分。对链内（folding）与链间（binding）分别进行上述次数的采样，并根据结果统计直方图信息。    
      3       1        0.000              ; dbin_f, dbin_b, vbin0    
链内与链间的采样间隔。    
上述采样的反应坐标是接触数Q。由于该模型依据Q来定义蛋白质的状态，因此这项采样至关重要。下述几项结果的采样并非必定进行，根据实际需要选择是否采用。    
          0         100         7.0   -350.0    ; IsEbin, nEbin, dEbin, vEbin0    
总能量采样。通过第一个参数0或1判断是否对该项结果进行采样。后三个参数分别为采样次数、采样间隔能量值与起始能量值。    
          1          100         0.6   -57    ; IsEbbin, nEbbin, dEbbin, vEbbin0    
结合能采样。定义同上。    
          0          100        1.000     0.000    ; IsRbin, nRbin, dRbin, vRbin0    
回转半径采样。定义同上。    
          1          100        0.5   -25  0.001   ; IsWbin, nWbin, dWbin, vWbin0, cri_Qb    
Qw采样。定义同上。Qw为改良过的结合界面接触数Qb。    
      7500.0      10.000                           ; PBC: pL, dl    
周期性边界条件。    
      1.200       1.200         0.96              ; Alpha1, Alpha2, Beta    
相互作用强度参数，分别定义两条链内和链间的相互作用强度。该模型把溶剂环境改变对蛋白质相互作用的影响用上述三个参数来描述，调整该参数即可模拟不同的溶剂条件。    
      0.000       5.000                           ; Delta, CritR_non    
非天然相互作用判定参数。    
