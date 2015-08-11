elfe_msgp.F90:343:    write(*,'(i4,2a)') myrank,': ABORT: ',string
elfe_msgp.F90:344:    if(lopen) write(11,'(i4,2a)') myrank,': ABORT: ',string
elfe_msgp.F90:349:      write(*,'(i4,2a)') myrank,': MPI ERROR: ',s
elfe_msgp.F90:350:      if(lopen) write(11,'(i4,2a)') myrank,': MPI ERROR: ',s
elfe_msgp.F90:406:      write(errmsg,*) 'PARALLEL_RRSYNC: recv start msg: ',idrcv,itrcv
elfe_msgp.F90:421:      write(errmsg,*) 'PARALLEL_RRSYNC: send start msg: ',idsnd,itsnd
elfe_msgp.F90:433:      write(errmsg,*) 'PARALLEL_RRSYNC: recv start msg: ',idrcv,itrcv
elfe_msgp.F90:445:        write(errmsg,*) 'PARALLEL_RRSYNC: send start msg: ',idsnd,itsnd
elfe_msgp.F90:482:  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
elfe_msgp.F90:557:  write(10,'(a,i8)') 'Number of neighbors:',nnbr
elfe_msgp.F90:558:  write(10,'(a)') '##########################################################'
elfe_msgp.F90:559:  write(10,'(a)') 'Element Receive Table:'
elfe_msgp.F90:561:    write(10,'(a,3i8)') 'nbrindx,rank,nerecv: ',&
elfe_msgp.F90:564:      write(10,'(t1,2i8)') ierecv(j,i),iegrecv(j,i)
elfe_msgp.F90:567:  write(10,'(a)') '##########################################################'
elfe_msgp.F90:646:  write(10,'(a)') 'Element Send Table:'
elfe_msgp.F90:648:    write(10,'(a,3i8)') 'nbrindx,rank,nesend: ',&
elfe_msgp.F90:651:      write(10,'(t1,2i8)') iesend(j,i),iegsend(j,i)
elfe_msgp.F90:654:  write(10,'(a)') '##########################################################'
elfe_msgp.F90:719:      write(errmsg,*) 'comm_table: error1 in ipgl: ',myrank,ip,iplg(ip)
elfe_msgp.F90:809:  write(10,'(a)') 'Node Receive Table:'
elfe_msgp.F90:811:    write(10,'(a,3i8)') 'nbrindx,rank,nprecv: ',&
elfe_msgp.F90:814:      write(10,*) 'Zero recv'
elfe_msgp.F90:815:!      write(errmsg,*) 'MSGP: Zero recv; see ctb*'
elfe_msgp.F90:819:      write(10,'(t1,2i8)') iprecv(j,i),ipgrecv(j,i)
elfe_msgp.F90:822:  write(10,'(a)') '##########################################################'
elfe_msgp.F90:896:  write(10,'(a)') 'Node Send Table:'
elfe_msgp.F90:898:    write(10,'(a,3i8)') 'nbrindx,rank,npsend: ',&
elfe_msgp.F90:901:      write(10,*) 'Zero send'
elfe_msgp.F90:902:!      write(errmsg,*) 'MSGP: Zero send; see ctb*'
elfe_msgp.F90:906:      write(10,'(t1,2i8)') ipsend(j,i),ipgsend(j,i)
elfe_msgp.F90:909:  write(10,'(a)') '##########################################################'
elfe_msgp.F90:919:            write(errmsg,*)'MSGP: address clash:',i,j,k,l,myrank
elfe_msgp.F90:953:      write(errmsg,*) 'comm_table: error1 in isgl: ',myrank,isd,islg(isd)
elfe_msgp.F90:1042:  write(10,'(a)') 'Side Receive Table:'
elfe_msgp.F90:1044:    write(10,'(a,3i8)') 'nbrindx,rank,nsrecv: ',&
elfe_msgp.F90:1047:      write(10,*)'Zero recv side'
elfe_msgp.F90:1048:!      write(errmsg,*) 'MSGP: Zero recv side; see ctb*'
elfe_msgp.F90:1052:      write(10,'(t1,4i8)') isrecv(j,i),isgrecv(j,i),iplg(isidenode(1:2,isrecv(j,i)))
elfe_msgp.F90:1055:  write(10,'(a)') '##########################################################'
elfe_msgp.F90:1121:  write(10,'(a)') 'Side Send Table:'
elfe_msgp.F90:1123:    write(10,'(a,3i8)') 'nbrindx,rank,nssend: ',&
elfe_msgp.F90:1126:      write(10,*)'Zero send side'
elfe_msgp.F90:1127:!      write(errmsg,*) 'MSGP: Zero send side; see ctb*'
elfe_msgp.F90:1131:      write(10,'(t1,4i8)') issend(j,i),isgsend(j,i),iplg(isidenode(1:2,issend(j,i)))
elfe_msgp.F90:1134:  write(10,'(a)') '##########################################################'
elfe_msgp.F90:1144:            write(errmsg,*)'MSGP: address clash (2):',i,j,k,l,myrank
grid_subs.F90:29:!    write(errmsg,*)'LINDEX: ',node,' is not in element ',ie
grid_subs.F90:105:  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
grid_subs.F90:107:  write(10,'(a,4i10)') '#',nea,ne,neg
grid_subs.F90:121:      write(10,'(a,2i8,4e14.6)') 'Element ',ie,ielg(ie),xctr(ie),yctr(ie),zctr(ie),dpe(ie)
grid_subs.F90:123:      write(10,'(a,2i8,4e14.6)') '# Element ',ie,ielg(ie),xctr(ie),yctr(ie),zctr(ie),dpe(ie)
grid_subs.F90:125:    write(10,'(a,3i8)') '####NODE:  ',(iplg(elnode(k,ie)),k=1,3)
grid_subs.F90:131:          write(errmsg,*)'Resident element having wrong nbr:',ie,ielg(ie),myrank
grid_subs.F90:141:        write(errmsg,*)'Check elnode or elside:',ielg(ie),(elnode(ip,ie),elside(ip,ie),ip=1,3)
grid_subs.F90:146:    write(10,'(a,3i8)') '####IC3:   ',(ibuf3(k),k=1,3)
grid_subs.F90:147:    write(10,'(a,3i8)') '####JS:    ',(islg(elside(k,ie)),k=1,3)
grid_subs.F90:148:    write(10,'(a,3i8)') '####SSIGN: ',(int(ssign(k,ie)),k=1,3)
grid_subs.F90:149:    write(10,'(a,64i8)') '####PList:',(ibuf1(k),ibuf2(k),k=1,j)
grid_subs.F90:156:  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
grid_subs.F90:158:  write(10,'(a,4i10)') '#',npa,np,npg
grid_subs.F90:178:          write(errmsg,*)'Surrounding element outside:',indel(k,ip),iplg(ip),k
grid_subs.F90:186:        write(errmsg,*)'Surrounding element not exist:',indel(k,ip),iplg(ip),k
grid_subs.F90:191:      write(10,'(a,2i8,4e14.6,2i4,50(i8,i4))') 'Node ',ip,iplg(ip),xnd(ip),ynd(ip),znd(ip),dp(ip), &
grid_subs.F90:194:      write(10,'(a,2i8,4e14.6,2i4,50(i8,i4))') '# Node ',ip,iplg(ip),xnd(ip),ynd(ip),znd(ip),dp(ip), &
grid_subs.F90:197:    write(10,'(a,64i8)') '####PList:',(ibuf1(k),ibuf2(k),k=1,j)
grid_subs.F90:204:          write(errmsg,*)'Surrounding node outside:',indnd(k,ip),iplg(ip),k
grid_subs.F90:210:        write(errmsg,*)'Surrounding node not exist:',indnd(k,ip),iplg(ip),k
grid_subs.F90:214:    write(10,'(a,64i8)')'Nbr nodes:',(ibuf3(k),k=1,nnp(ip))
grid_subs.F90:222:  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
grid_subs.F90:224:  write(10,'(a,4i10)') '#',nsa,ns,nsg
grid_subs.F90:233:      write(errmsg,*)'isdel(1,:) =0:',ngb1,ngb2
grid_subs.F90:247:      write(10,'(a,6i8,7e14.6,i4)') 'Side ',isd,isdgb,ngb1,ngb2,iegb1,iegb2, &
grid_subs.F90:250:      write(10,'(a,6i8,7e14.6,i4)') '# Side', isd,isdgb,ngb1,ngb2,iegb1,iegb2, &
grid_subs.F90:253:    write(10,'(a,64i8)') '####PList:',(ibuf1(k),ibuf2(k),k=1,j)
grid_subs.F90:260:  write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
grid_subs.F90:262:  write(10,'(a,i10)') 'Open bnd:',nope
grid_subs.F90:264:    write(10,*)'open bnd #',i,iopelg(i),(iplg(iond(i,j)),j=1,nond(i))
grid_subs.F90:266:  write(10,'(a,i10)') 'Land bnd:',nland
grid_subs.F90:268:    write(10,*)'land bnd #',i,(iplg(ilnd(i,j)),j=1,nlnd(i))
grid_subs.F90:277:    write(32,'(i8,1x,i4)')(ie,iegrpv(ie),ie=1,ne_global)
grid_subs.F90:371:        write(errmsg,*)'PARTITION: remove south pole:',i,rearth+znd(i)
grid_subs.F90:576:!  if(myrank==0) write(*,'(/a)') 'ParMETIS Partitioning:'
grid_subs.F90:669:        write(errmsg,*)'Wrong kz:',kz
grid_subs.F90:673:        write(errmsg,*)'h_s needs to be larger:',h_s
grid_subs.F90:686:          write(errmsg,*)'Illegal Z level:',k
grid_subs.F90:690:          write(errmsg,*)'z-level inverted:',k
grid_subs.F90:702:        write(errmsg,*)'h_c needs to be larger:',h_c
grid_subs.F90:706:        write(errmsg,*)'Wrong theta_b:',theta_b
grid_subs.F90:710:        write(errmsg,*)'Wrong theta_f:',theta_f
grid_subs.F90:723:          write(errmsg,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
grid_subs.F90:755:!    write(10,*)'Sample z coordinates'
grid_subs.F90:757:!    write(10,*)'h_c= ',h_c,', h_s=',h_s
grid_subs.F90:759:!      write(10,*)'Depth= ',buf1(i)
grid_subs.F90:768:!        write(10,*)k,zz
grid_subs.F90:865:      write(errmsg,*) 'AQUIRE_HGRID: Unknown type of element',iegb,j
grid_subs.F90:887:        write(11,*) 'Hanging node:',ipgb
grid_subs.F90:899:  if(myrank==0.and..not.full_aquire) write(16,*)'mnei= ',mnei
grid_subs.F90:933:            write(errmsg,*) 'Triangles ', iegb, ' and ', jegb, ' have opposite orientation'
grid_subs.F90:961:            write(errmsg,'(a,10i6)') 'AQUIRE_HGRID: Wrong ball info', &
grid_subs.F90:971:    write(errmsg,*) &
grid_subs.F90:978:      write(16,'(/a,4i10)')'Global Grid Size (ne,np,ns,nvrt): ',ne_global,np_global,ns_global,nvrt
grid_subs.F90:1011:    write(errmsg,*) 'neta_global /= total # of open bnd nodes',neta_global,nt
grid_subs.F90:1039:      write(errmsg,*) 'Looped open bnd:',k
grid_subs.F90:1067:    write(errmsg,*) 'AQUIRE_HGRID: nvel_global /= total # of land bnd nodes', &
grid_subs.F90:1128:        write(errmsg,*)'Illegal bnd node',i
grid_subs.F90:1150:        write(errmsg,*)'Incomplete ball',i
grid_subs.F90:1159:!        write(errmsg,*)'Failed to find local index:',i,new
grid_subs.F90:1167:!          write(errmsg,*)'Broken ball:',i
grid_subs.F90:1173:!          write(errmsg,*)'Too many neighbor nodes:',i,mnei
grid_subs.F90:1348:            write(errmsg,*)'Overflow in ieg:',npi*mnei
grid_subs.F90:1381:              write(errmsg,*)'Overflow in ipg:',neg*3
grid_subs.F90:1401:              write(errmsg,*)'Overflow in isg:',neg*3
grid_subs.F90:1441:      write(16,'(/a)') '**********Augmented Subdomain Sizes**********'
grid_subs.F90:1442:      write(16,'(10a)') ' rank', &
grid_subs.F90:1447:        write(16,'(i5,9i8)') i, &
grid_subs.F90:1594:        write(errmsg,*)'Too many neighbor nodes (0):',i,mnei+1
grid_subs.F90:1605:            write(errmsg,*)'Broken ball:',i
grid_subs.F90:1611:            write(errmsg,*)'Too many neighbor nodes:',i,mnei+1
grid_subs.F90:1756:        write(errmsg,*)'AQUIRE_HGRID: orientation wrong:',ielg(ie),egb1,&
grid_subs.F90:1768:        write(errmsg,*)'AQUIRE_HGRID: axes wrong',ielg(ie),xtmp,ytmp,dptmp
grid_subs.F90:1785:      write(errmsg,'(a,2i8)') 'AQUIRE_HGRID: negative area at',ie,ielg(ie)
grid_subs.F90:1794:      write(16,*)'Max. dot product of 3 axes=',real(dptmp) !thetan
grid_subs.F90:1841:        write(errmsg,*) 'AQUIRE_HGRID: Zero side',jsj
grid_subs.F90:1857:        write(errmsg,*)'First element empty:',ie,iegb,j
grid_subs.F90:1897:        write(errmsg,*)'AQUIRE_HGRID: 0 ys-vector',iplg(isidenode(1:2,j))
grid_subs.F90:1907:        write(errmsg,*)'AQUIRE_HGRID: 0 zs-vector',iplg(isidenode(1:2,j))
grid_subs.F90:1923:        write(errmsg,*)'AQUIRE_HGRID: 0 xs-vector',iplg(isidenode(1:2,j))
grid_subs.F90:1941:      !  write(12,*)'sample sframe:',iplg(n1),iplg(n2),xtmp,ytmp,sframe(:,:,j)
grid_subs.F90:1949:      write(16,*)'Max. deviation between ze and zs axes=',real(xtmp) !thetan
grid_subs.F90:1950:      write(16,*)'Max. dot prod. between ys and zs axes=',real(ytmp) !realvalue
grid_subs.F90:2027:        write(errmsg,*)'AQUIRE_HGRID: nond(nope)>mnond',nond(nope),mnond
grid_subs.F90:2145:          write(errmsg,*)'agquire_hgrid: node on more than 2 open bnds:',ipgb
grid_subs.F90:2189:              write(errmsg,*)'aquire_hgrid: impossible (1)',n1,n2,isdel(1:2,k),ielg(isdel(1:2,k)),ielg(ie),iplg(isidenode(1:2,isd))
grid_subs.F90:2209:    write(16,'(/a)') '**********Global Boundary Sizes**********'
grid_subs.F90:2210:    write(16,'(4a)') '    nope','    neta','   nland','    nvel'
grid_subs.F90:2211:    write(16,'(4i8)') nope_global,neta_global,nland_global,nvel_global
grid_subs.F90:2221:    write(16,'(/a)') '**********Augmented Subdomain Boundary Sizes**********'
grid_subs.F90:2222:    write(16,'(5a)') '    rank','    nope','    neta','   nland','    nvel'
grid_subs.F90:2224:      write(16,'(5i8)') i,irbuf(4*i+1),irbuf(4*i+2),irbuf(4*i+3),irbuf(4*i+4)
grid_subs.F90:2226:    write(16,*)
io_subs.F90:25:      write(32,*) 'Sidegrid'
io_subs.F90:26:      write(32,*) ns
io_subs.F90:29:          write(32,*) i,real(xcj(i)),real(ycj(i)),real(dps(i))
io_subs.F90:33:          write(32,*) i,real(tmp1),real(tmp2),real(dps(i))
io_subs.F90:40:      write(32,*) 'centers pts'
io_subs.F90:41:      write(32,*) ne
io_subs.F90:44:          write(32,*) i,real(xctr(i)),real(yctr(i)),real(dpe(i))
io_subs.F90:48:          write(32,*) i,real(tmp1),real(tmp2),real(dpe(i))
io_subs.F90:110:        write(36,'(2a)') '# ','********** Timer Index Mapping **********'
io_subs.F90:111:        write(36,'(2a)') '# ','00 -- Total'
io_subs.F90:112:        write(36,'(2a)') '# ','01 -- Init Section'
io_subs.F90:113:        write(36,'(2a)') '# ','02 -- Timestepping Section'
io_subs.F90:114:        write(36,'(2a)') '# ','03 -- Forcings & Prep Section'
io_subs.F90:115:        write(36,'(2a)') '# ','04 -- Backtracking Section'
io_subs.F90:116:        write(36,'(2a)') '# ','05 -- Turbulence Closure Section'
io_subs.F90:117:        write(36,'(2a)') '# ','06 -- Matrix Preparation Section'
io_subs.F90:118:        write(36,'(2a)') '# ','07 -- Wave-Cont. Solver Section'
io_subs.F90:119:        write(36,'(2a)') '# ','08 -- Momentum Eqs. Solve Section'
io_subs.F90:120:        write(36,'(2a)') '# ','09 -- Transport Section'
io_subs.F90:121:        write(36,'(2a)') '# ','10 -- Recomputing Levels Section'
io_subs.F90:122:        write(36,'(2a)') '# ','11 -- Conservation Check Section'
io_subs.F90:123:        write(36,'(2a)') '# ','12 -- Global Output Section'
io_subs.F90:124:        write(36,'(2a)') '# ','13 -- Hotstart Section'
io_subs.F90:128:        write(36,'(/)')
io_subs.F90:129:        write(36,'(2a)') '# ','********** Average, Min & Max Times in secs **********'
io_subs.F90:131:        write(36,'(11a)') 'ID', &
io_subs.F90:135:          write(36,'(i2.2,2(e13.6,2(e13.6,i13)))') i, &
io_subs.F90:154:        write(36,'(/)')
io_subs.F90:155:        write(36,'(a)') '# ********** Computation Times (sec) For Each MPI Task **********'
io_subs.F90:156:        write(36,'(a)') '# ********** Rows = Ranks; Columns = Timers      **********'
io_subs.F90:157:        write(36,'(a,20i13)') '# Rank',(i,i=0,13)
io_subs.F90:161:      write(36,'(a,i4.4,20e13.6)') '# ',myrank,(wtimer(i,1),i=0,13)
io_subs.F90:175:        write(36,'(/)')
io_subs.F90:176:        write(36,'(a)') '# ********** Communication Times For Each MPI Task **********'
io_subs.F90:177:        write(36,'(a)') '# ********** Rows = Ranks; Columns = Timers        **********'
io_subs.F90:178:        write(36,'(a,20i13)') '# Rank',(i,i=0,13)
io_subs.F90:182:      write(36,'(a,i4.4,20e13.6)') '# ',myrank,(wtimer(i,2),i=0,13)
io_subs.F90:268:            if(myrank==0) write(86,*)varname,' = ',varvalue2
io_subs.F90:273:            if(myrank==0) write(86,*)varname,' = ',ivarvalue
io_subs.F90:278:            if(myrank==0) write(86,*)varname,' = ',real(varvalue1)
io_subs.F90:281:            write(errmsg,*)'read_param: unknown type:',vartype
io_subs.F90:293:      write(errmsg,*)'Failed to find parameter:',varname
misc_subs.F90:58:        write(errmsg,*)'ZCOOR: dry location:',dp(inode),eta2(inode),itag
misc_subs.F90:74:            write(errmsg,*)'ZCOOR: Pls choose a larger h_c:',eta2(inode),h_c,itag
misc_subs.F90:94:            write(errmsg,*)'ZCOOR: Cannot find a bottom level:',dp(inode),itag
misc_subs.F90:105:!          write(errmsg,*)'ZCOOR: elev<hsm:',eta,itag
misc_subs.F90:123:          write(12,*)'ZCOOR: Inverted z-level:',itag,ivcor,k,kbpl,iplg(inode),eta2(inode),dp(inode),ztmp(k),ztmp(k-1),sigma_lcl(kbpl:nvrt,inode)
misc_subs.F90:124:          write(errmsg,*)'ZCOOR: Inverted z-level:',itag,ivcor,k,kbpl,iplg(inode),eta2(inode),dp(inode),ztmp(k),ztmp(k-1)
misc_subs.F90:201:!        write(12,*)'it=',it  
misc_subs.F90:205:!          write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
misc_subs.F90:207:!          write(10,*)np
misc_subs.F90:209:!            write(10,*)iplg(i),real(eta2(i))
misc_subs.F90:211:!          write(10,*)ns
misc_subs.F90:213:!            write(10,*)i,iplg(isidenode(1:2,i)),real(su2(nvrt,i)),real(sv2(nvrt,i))
misc_subs.F90:253:            if(myrank==0) write(16,*)'doing final extrapolation in levels1...'
misc_subs.F90:261:                write(errmsg,*)'LEVELS1: bnd side (2):',isdel(:,isd),iplg(isidenode(1:2,isd))
misc_subs.F90:319:                        write(errmsg,*)'LEVELS1: Failed to wet element:',ielg(ie2),iplg(nodeA)
misc_subs.F90:376:!                write(12,*)'Make dry:',itr,iplg(nd)
misc_subs.F90:395:              write(errmsg,*)'LEVELS1: bnd side:',isdel(:,isd),iplg(isidenode(1:2,isd))
misc_subs.F90:412:              write(errmsg,*)'Frontier node outside, or on the interface:', &
misc_subs.F90:415:              write(12,*)'LEVELS1: fatal error message'
misc_subs.F90:418:                  write(12,*)l,iplg(isidenode(1:2,l))
misc_subs.F90:419:                  write(12,*)l,ielg(isdel(1:2,l)),idry_e2(isdel(1:2,l)),idry_e(isdel(1:2,l))
misc_subs.F90:423:                write(12,*)l,idry_e2(l),idry_e(l)
misc_subs.F90:433:                  write(errmsg,*)'Failed to wet element (13):',ielg(ie),iplg(nd),iplg(nodeA)
misc_subs.F90:439:!              write(12,*)'Make wet:',itr,iplg(nodeA),ielg(ie)
misc_subs.F90:602:          write(16,*)'see fort.7 for # of iterations used in LEVELS1...'
misc_subs.F90:603:          write(7,*)it,itr
misc_subs.F90:652:          write(errmsg,*)'Deep depth dry:',iplg(i)
misc_subs.F90:662:          !  write(errmsg,*)'levels1: (2):',i,dp(i)+eta2(i)
misc_subs.F90:674:!              write(errmsg,*)'Pls choose a larger h_c (2):',eta2(i),h_c
misc_subs.F90:694:!                write(errmsg,*)'Cannot find a bottom level for node (3):',i
misc_subs.F90:701:!              write(errmsg,*)'Impossible 92:',kbp(i),kz,i
misc_subs.F90:712:!              write(errmsg,*)'Inverted z-levels at:',i,k,znl(k,i)-znl(k-1,i),eta2(i),hmod(i)
misc_subs.F90:722:!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
misc_subs.F90:725:!      write(10,*)'Time step=',it
misc_subs.F90:726:!      write(10,*)'Node'
misc_subs.F90:728:!        write(10,*)i,iplg(i),dp(i),eta2(i)
misc_subs.F90:739:          write(errmsg,*)'level1: Element-node inconsistency (0):',ielg(i),idry_e(i), &
misc_subs.F90:747:            write(errmsg,*)'Weird element:',k,i,ze(k,i),ze(k-1,i)
misc_subs.F90:761:            write(errmsg,*)'Side-node inconsistency:',it,islg(i),'node:',iplg(n1),iplg(n2), &
misc_subs.F90:767:            write(errmsg,*)'Weird side:',islg(i),iplg(n1),iplg(n2),eta2(n1),eta2(n2)
misc_subs.F90:774:              write(errmsg,*)'Weird side:',k,iplg(n1),iplg(n2),znl(max(k,kbp(n1)),n1), &
misc_subs.F90:815:!                  write(12,*)'Isolated rewetted node:',iplg(i)
misc_subs.F90:867:!                  write(12,*)'Isolated rewetted side:',i,iplg(n1),iplg(n2)
misc_subs.F90:997:            write(errmsg,*)'Deep depth dry:',i
misc_subs.F90:1013:!              write(errmsg,*)'Pls choose a larger h_c (1):', ' node:', iplg(i), ', elev prev:', eta1(i), ', elev cur:', eta2(i), ', h_c:', h_c
misc_subs.F90:1033:!                write(errmsg,*)'Cannot find a bottom level for node (3):',i
misc_subs.F90:1040:!              write(errmsg,*)'Impossible 92:',kbp(i),kz,i
misc_subs.F90:1051:!              write(errmsg,*)'Inverted z-levels at:',i,k,znl(k,i)-znl(k-1,i),eta2(i),hmod(i)
misc_subs.F90:1061:!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
misc_subs.F90:1064:!      write(10,*)'Time step=',it
misc_subs.F90:1065:!      write(10,*)'Node'
misc_subs.F90:1067:!        write(10,*)i,iplg(i),dp(i),eta2(i),idry2(i)
misc_subs.F90:1080:!      write(10,*)'Element'
misc_subs.F90:1082:!        write(10,*)i,ielg(i),idry_e2(i)
misc_subs.F90:1104:!      write(10,*)'nodes'
misc_subs.F90:1106:!        write(10,*)i,iplg(i),idry2(i),np
misc_subs.F90:1115:!          write(errmsg,*)'levels0: weird wet node:',iplg(i),eta2(i),dp(i),idry2(i)
misc_subs.F90:1128:!          write(errmsg,*)'Node-element inconsistency:',iplg(i),idry2(i),(idry_e2(indel(j,i)),j=1,nne(i))
misc_subs.F90:1142:          write(errmsg,*)'level0: Element-node inconsistency (0):',ielg(i),idry_e2(i), &
misc_subs.F90:1150:            write(errmsg,*)'Weird element:',k,i,ze(k,i),ze(k-1,i)
misc_subs.F90:1187:!                  write(12,*)'Isolated rewetted node:',iplg(i)
misc_subs.F90:1220:!      write(10,*)'Side'
misc_subs.F90:1222:!        write(10,*)i,islg(i),idry_s2(i),ns
misc_subs.F90:1233:!            write(errmsg,*)'Element-side inconsistency:',ielg(i),islg(isd),idry_s2(isd)
misc_subs.F90:1250:!          write(errmsg,*)'Side-element inconsistency:',islg(i),idry_s2(i), &
misc_subs.F90:1264:            write(errmsg,*)'Side-node inconsistency:',it,islg(i),'node:',iplg(n1),iplg(n2), &
misc_subs.F90:1270:            write(errmsg,*)'Weird side:',islg(i),iplg(n1),iplg(n2),eta2(n1),eta2(n2)
misc_subs.F90:1277:              write(errmsg,*)'Weird side:',k,iplg(n1),iplg(n2),znl(max(k,kbp(n1)),n1), &
misc_subs.F90:1336:!                  write(12,*)'Isolated rewetted side:',i,iplg(n1),iplg(n2)
misc_subs.F90:1532:            write(errmsg,*)'Isolated wet node (8):',i
misc_subs.F90:1647:            write(errmsg,*)'nodalvel: Isolated open bnd node:',iplg(i),isbnd(1:2,i)
misc_subs.F90:1697:!	      write(*,*)'Wrong element ball'
misc_subs.F90:1745:        write(errmsg,*)'k1>k2 in vinter()'
misc_subs.F90:1776:            write(errmsg,*)'Failed to find a level in vinter():',kout,zt,(za(k),k=k1,k2)
misc_subs.F90:1831:        write(errmsg,*)'EQSTATE: Impossible dry (7):',tem,sal
misc_subs.F90:1837:          write(12,*)'Invalid temp. or salinity for density:',tem,sal
misc_subs.F90:1851:        write(errmsg,*)'Weird density:',eqstate,tem,sal
misc_subs.F90:1857:!      if (myrank==0) write(16,*)'sediment density effect'
misc_subs.F90:1861:!        write(12,*)Srho(ised),sconc(ised),eqstate
misc_subs.F90:1865:!        write(12,*)SedDen,eqstate
misc_subs.F90:1875:            write(errmsg,*)'EQSTATE: Impossible (8):',indx,ised,rho_w,Srho(ised),sconc(ised)
misc_subs.F90:1900:        write(errmsg,*)'Wrong input level:',j
misc_subs.F90:1932:        write(errmsg,*)'Unknown ASM:',mid
misc_subs.F90:1977:        write(errmsg,*)'Check inputs in rint_lag:',Nmin,Nmax
misc_subs.F90:1981:        write(errmsg,*)'Wrong k:',k
misc_subs.F90:1985:        write(errmsg,*)'m<1',m
misc_subs.F90:1989:        write(errmsg,*)'m>3 not covered presently' 
misc_subs.F90:1993:        write(errmsg,*)'Re-dimension sigmap'
misc_subs.F90:2001:         write(errmsg,*)'Weird indices:',j1,j2
misc_subs.F90:2018:          write(errmsg,*)'Miscount:',id,j2-j1,m
misc_subs.F90:2074:          write(errmsg,*)'Not covered:',id
misc_subs.F90:2084:          write(errmsg,*)'sigma_prod index out of bound (2)'
misc_subs.F90:2127:        write(errmsg,*)'Negative hh in covar:',hh
misc_subs.F90:2145:        write(errmsg,*)'Unknown covariance function option:',kr_co
misc_subs.F90:2185:!        write(errmsg,*)'EVAL_CUBIC: xmin>xmax:',xmin,xmax
misc_subs.F90:2213:          write(errmsg,*)'EVAL_CUBIC: Falied to find:',i,xtmp,xmin,xmax
misc_subs.F90:2247:            write(errmsg,*)'CUBIC_SP: bottom problem:',xcor(k+1),xcor(k)
misc_subs.F90:2255:            write(errmsg,*)'CUBIC_SP: surface problem:',xcor(k),xcor(k-1)
misc_subs.F90:2265:            write(errmsg,*)'CUBIC_SP: middle problem:',xcor(k),xcor(k-1),xcor(k+1)
misc_subs.F90:2538:            write(errmsg,*)'hgrad_nodes: node3 dry',iplg(node3),ielg(ie)
misc_subs.F90:2555:              write(errmsg,*)'hgrad_nodes: no index found:',iplg(nd),ielg(ie)
misc_subs.F90:2662:            write(errmsg,*)'hgrad_nodes failure:',iplg(node1),iplg(node2)
misc_subs.F90:2780:        write(errmsg,*)'COMPUTE_LL: rad=0:',xg,yg,zg,rad
misc_subs.F90:2931:        write(errmsg,*)'WBL: check inputs:',z0,ubm,wfr
misc_subs.F90:2961:        write(errmsg,*)'WBL: exponent too large (1):',tmp,rkn,wfr,c_mu,ubm
misc_subs.F90:2975:!            write(*,*)'wave bottom layer did not converge:',rmu,rmu2,tau_wm,fw,ubm,phi_cw,c_mu,wfr
misc_subs.F90:2984:            write(errmsg,*)'WBL: exponent too large (2):',tmp,rkn,wfr,c_mu,ubm,rmu2,rmu
misc_subs.F90:3067:          write(errmsg,*)'AREA_COORD: failed to fix',arco(1:3)
solver_subs.F90:105:    write(errmsg,*)'JCG: 0 initial error:',rdotr0
solver_subs.F90:111:    write(33,'(//a,i8)') '********CG Solve at timestep ',itime
solver_subs.F90:112:    write(33,'(a,i6,2e14.6)') &
solver_subs.F90:119:      if(myrank==0) write(33,*)'JCG converged in ',itn,' iterations'
solver_subs.F90:123:      if(myrank==0) write(33,*)'JCG did not converge in ',mxitn,' iterations' 
solver_subs.F90:211:    if(mod(itn,moitn)==0.and.myrank==0) write(33,'(a,i6,2e14.6)') &
solver_subs.F90:303:          write(errmsg,*)'Not symmetric:',iplg(ip),k,l,blockj(k,l,ip)-blockj(l,k,ip)
solver_subs.F90:307:      !write(12,*)1/qmatr(k,0,0,ip),ip,k,blockj(k,kbp_e(ip):(nvrt-1),ip) 
solver_subs.F90:316:    write(29,'(//a,i8)') '********CG2 Solve at timestep ',itime
solver_subs.F90:317:    write(29,*)'done pre-conditioner'
solver_subs.F90:318:    if(rat_max2_gb>0) write(16,*)'Max. vertical/horizontal ratio and threshold=',sqrt(rat_max2_gb),threshold_rat
solver_subs.F90:358:    write(29,'(a,i6,3e14.6)')'Itn, 2Norm, Rnorm, tol2: ',itn,rdotr,rdotr/rdotr0,rtol2
solver_subs.F90:364:      if(myrank==0) write(29,*)'JCG2 converged in ',itn,' iterations'
solver_subs.F90:368:      if(myrank==0) write(29,*)'JCG2 did not converge in ',mxitn,' iterations'
solver_subs.F90:480:    if(mod(itn,moitn)==0.and.myrank==0) write(29,'(a,i6,3e14.6)') &
wwm_airsea.F90:104:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111115,'(2I10,5F15.8)') I, ITAUMAX, XI, TAUW(IJ), DELTAUW 
wwm_airsea.F90:116:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111115,'(5F15.8)') UST2, ARG, TAUW(IJ), EPS1,  Z0(IJ)
wwm_ardhuin_new.F90:228:!        WRITE(5001,*) 'DTH'
wwm_ardhuin_new.F90:229:!        WRITE(5001,*) DTH
wwm_ardhuin_new.F90:230:!        WRITE(5001,*) 'FR1'
wwm_ardhuin_new.F90:231:!        WRITE(5001,*) FR1
wwm_ardhuin_new.F90:232:!        WRITE(5001,*) 'TH'
wwm_ardhuin_new.F90:233:!        WRITE(5001,*) TH
wwm_ardhuin_new.F90:234:!        WRITE(5001,*) 'ESIN, ECOS, EC2'
wwm_ardhuin_new.F90:235:!        WRITE(5001,*) ESIN, ECOS, EC2
wwm_ardhuin_new.F90:239:!        WRITE(5001,*) 'WNMEANP, WNMEANPTAIL'
wwm_ardhuin_new.F90:240:!        WRITE(5001,*) WNMEANP, WNMEANPTAIL
wwm_ardhuin_new.F90:247:!        WRITE(5001,*) 'XFR, SIGMA, SXFR'
wwm_ardhuin_new.F90:248:!        WRITE(5001,*) XFR, SIGMA, SXFR
wwm_ardhuin_new.F90:256:!        WRITE(5001,*) 'SIGMA'
wwm_ardhuin_new.F90:257:!        WRITE(5001,*)  SIGMA
wwm_ardhuin_new.F90:258:!        WRITE(5001,*) 'SIG'
wwm_ardhuin_new.F90:259:!        WRITE(5001,*)  SIG
wwm_ardhuin_new.F90:260:!        WRITE(5001,*) 'DSIP'
wwm_ardhuin_new.F90:261:!        WRITE(5001,*)  DSIP
wwm_ardhuin_new.F90:277:!        WRITE(5001,*) 'SIG2'
wwm_ardhuin_new.F90:278:!        WRITE(5001,*) SIG2
wwm_ardhuin_new.F90:279:!        WRITE(5001,*) 'DSII'
wwm_ardhuin_new.F90:280:!        WRITE(5001,*) DSII
wwm_ardhuin_new.F90:281:!        WRITE(5001,*) 'DDEN'
wwm_ardhuin_new.F90:282:!        WRITE(5001,*) DDEN
wwm_ardhuin_new.F90:283:!        WRITE(5001,*) 'DDEN2'
wwm_ardhuin_new.F90:284:!        WRITE(5001,*) DDEN2
wwm_ardhuin_new.F90:296:!        WRITE(5001,*) 'FTE, FTF, FACHF, FACHFE'
wwm_ardhuin_new.F90:297:!        WRITE(5001,*) FTE, FTF, FACHF, FACHFE
wwm_ardhuin_new.F90:307:!        WRITE(5001,*) 'SSWELLF'
wwm_ardhuin_new.F90:308:!        WRITE(5001,*) SSWELLF
wwm_ardhuin_new.F90:437:          WRITE(DBG%FHNDL,*) 'MSC AND MDC READ FROM FILE AND SET IN WWMINPUT.NML ARE NOT EQUAL -STOP-'
wwm_ardhuin_new.F90:438:          WRITE(DBG%FHNDL,*) MSC_TEST, MSC
wwm_ardhuin_new.F90:439:          WRITE(DBG%FHNDL,*) MDC_TEST, MDC 
wwm_ardhuin_new.F90:572:!          WRITE(DBG%FHNDL,*) IK, EB(IK), IK, ITH, A(ITH,IK)
wwm_ardhuin_new.F90:847:        !WRITE(DBG%FHNDL,'(A10,I10,5F15.8)') 'TEST IND',IND, XI, AORB, Z0NOZ, ABMIN, DELAB
wwm_ardhuin_new.F90:874:      !WRITE(DBG%FHNDL,*) 'TAU USTAR', TAUX, TAUY, UST, USDIR, USTAR
wwm_ardhuin_new.F90:896:        !WRITE(DBG%FHNDL,*) 'USTP', IK, USTP, STRESSSTAB(ISTAB,1), STRESSSTAB(ISTAB,2), TTAUWSHELTER
wwm_ardhuin_new.F90:926:        !WRITE(DBG%FHNDL,*) 'UCN', IK, IS, USTP, CM ,K(IK), DDEN2(IS), Z0
wwm_ardhuin_new.F90:939:            !WRITE(DBG%FHNDL,*) 'ZLOG', IK, ITH, ZCN, ZARG, X, KAPPA, UCN
wwm_ardhuin_new.F90:953:              !WRITE(DBG%FHNDL,*) 'DSTAB', DSTAB(ISTAB,IS), CONST,EXP(ZLOG),ZLOG**4,UCN**2,COSWIND,SSINTHP
wwm_ardhuin_new.F90:1024:         !WRITE(DBG%FHNDL,*) ITH, IS, A(IS), (MAX(COSWIND,ZERO))**3
wwm_ardhuin_new.F90:1030:      !WRITE(DBG%FHNDL,*) 'TAUPX', TAUPX, TAUPY, TTAUWSHELTER, XSTRESS, YSTRESS
wwm_ardhuin_new.F90:1173:!      WRITE(6,*) 'INSIN4:',FLTABS, SSDSDTH, SSDSC3, SSDSBCK
wwm_ardhuin_new.F90:1178:          WRITE(STAT%FHNDL,*) 'Computing 3D lookup table... please wait ...'
wwm_ardhuin_new.F90:1316:!        WRITE(6,*) 'INSIN4b:',DIKCUMUL                               
wwm_ardhuin_new.F90:1894:!      write(DBG%FHNDL,*) z0, ustar, windspeed
wwm_ardhuin_new.F90:2056:!WRITE(6,*) 'SDS IK:',IK,ASUM,MSSLONG(IK),BTH0(IK),FACSTRAIN
wwm_ardhuin_new.F90:2075:        !WRITE(DBG%FHNDL,*) 'FACSAT', FACSAT
wwm_ardhuin_new.F90:2076:        !WRITE(DBG%FHNDL,*) 'SATWEIGHTS', SATWEIGHTS 
wwm_ardhuin_new.F90:2077:        !WRITE(DBG%FHNDL,*) 'SATINDICES', SATINDICES
wwm_ardhuin_new.F90:2078:        !WRITE(DBG%FHNDL,*) 'SUMS', SUM(BTH), SUM(BTH0), SUM(A)
wwm_ardhuin_new.F90:2369:          !WRITE(DBG%FHNDL,*) 'TURBULENCE', IS, D(IS), (SSDSC(3)*RENEWALFREQ+DTURB)
wwm_ardhuin_old.F90:212:!        WRITE(5001,*) 'DTH'
wwm_ardhuin_old.F90:213:!        WRITE(5001,*) DTH
wwm_ardhuin_old.F90:214:!        WRITE(5001,*) 'FR1'
wwm_ardhuin_old.F90:215:!        WRITE(5001,*) FR1
wwm_ardhuin_old.F90:216:!        WRITE(5001,*) 'TH'
wwm_ardhuin_old.F90:217:!        WRITE(5001,*) TH
wwm_ardhuin_old.F90:218:!        WRITE(5001,*) 'ESIN, ECOS, EC2'
wwm_ardhuin_old.F90:219:!        WRITE(5001,*) ESIN, ECOS, EC2
wwm_ardhuin_old.F90:225:!        WRITE(5001,*) 'WNMEANP, WNMEANPTAIL'
wwm_ardhuin_old.F90:226:!        WRITE(5001,*) WNMEANP, WNMEANPTAIL
wwm_ardhuin_old.F90:233:!        WRITE(5001,*) 'XFR, SIGMA, SXFR'
wwm_ardhuin_old.F90:234:!        WRITE(5001,*) XFR, SIGMA, SXFR
wwm_ardhuin_old.F90:242:!        WRITE(5001,*) 'SIGMA'
wwm_ardhuin_old.F90:243:!        WRITE(5001,*)  SIGMA
wwm_ardhuin_old.F90:244:!        WRITE(5001,*) 'SIG'
wwm_ardhuin_old.F90:245:!        WRITE(5001,*)  SIG
wwm_ardhuin_old.F90:246:!        WRITE(5001,*) 'DSIP'
wwm_ardhuin_old.F90:247:!        WRITE(5001,*)  DSIP
wwm_ardhuin_old.F90:263:!        WRITE(5001,*) 'SIG2'
wwm_ardhuin_old.F90:264:!        WRITE(5001,*) SIG2
wwm_ardhuin_old.F90:265:!        WRITE(5001,*) 'DSII'
wwm_ardhuin_old.F90:266:!        WRITE(5001,*) DSII
wwm_ardhuin_old.F90:267:!        WRITE(5001,*) 'DDEN'
wwm_ardhuin_old.F90:268:!        WRITE(5001,*) DDEN
wwm_ardhuin_old.F90:269:!        WRITE(5001,*) 'DDEN2'
wwm_ardhuin_old.F90:270:!        WRITE(5001,*) DDEN2
wwm_ardhuin_old.F90:282:!        WRITE(5001,*) 'FTE, FTF, FACHF, FACHFE'
wwm_ardhuin_old.F90:283:!        WRITE(5001,*) FTE, FTF, FACHF, FACHFE
wwm_ardhuin_old.F90:293:!        WRITE(5001,*) 'SSWELLF'
wwm_ardhuin_old.F90:294:!        WRITE(5001,*) SSWELLF
wwm_ardhuin_old.F90:547:!          WRITE(DBG%FHNDL,*) IK, EB(IK), IK, ITH, A(ITH,IK)
wwm_ardhuin_old.F90:771:!/DEBUG    IF(HSBLOW.GT.1.00) WRITE(6,*) 'HS IN SIN3    :',IX,HSBLOW,IS,IK,ITH,A(5)
wwm_ardhuin_old.F90:798:          !WRITE(DBG%FHNDL,*) UORBT, AORB, EB, SIG(IK)**2, DDEN(IK), CG(IK)
wwm_ardhuin_old.F90:821:       !WRITE(DBG%FHNDL,*) 'SWELLFT', SWELLFT(IND),DELI2,SWELLFT(IND+1),DELI1, DELAB
wwm_ardhuin_old.F90:822:       !WRITE(DBG%FHNDL,*) 'URBOT', UORB, AORB, Z0NOZ,XI
wwm_ardhuin_old.F90:847:       !WRITE(DBG%FHNDL,*) 'TAU USTAR', TAUX, TAUY, UST, USDIR, USTAR
wwm_ardhuin_old.F90:867:        !WRITE(DBG%FHNDL,*) 'MESS', IK, UST, STRESSSTAB(ISTAB,1), STRESSSTAB(ISTAB,2), TTAUWSHELTER
wwm_ardhuin_old.F90:886:        !WRITE(DBG%FHNDL,*) 'UCN', IK, IS, USTP, CM ,K(IK), DDEN2(IS), Z0
wwm_ardhuin_old.F90:899:            !WRITE(DBG%FHNDL,*) 'ZLOG', IK, ITH, ZCN, ZARG, X, KAPPA, UCN
wwm_ardhuin_old.F90:912:              !WRITE(DBG%FHNDL,*) DSTAB(ISTAB,IS), CONST,EXP(ZLOG),ZLOG**4,UCN**2,COSWIND,SSINTHP
wwm_ardhuin_old.F90:967:        !WRITE(DBG%FHNDL,*) 'DSTAB', DSTAB(3,:)
wwm_ardhuin_old.F90:968:        !WRITE(DBG%FHNDL,*) 'STRESSTAB', STRESSSTAB (3,1), STRESSSTAB (3,2), STRESSSTABN(3,1), STRESSSTABN(3,2)
wwm_ardhuin_old.F90:969:        !WRITE(DBG%FHNDL,*) FW, UORB
wwm_ardhuin_old.F90:970:       !  WRITE(995,'(A,11G14.5)') 'NEGSTRESS:    ',TAUWNX,TAUWNY,FW*UORB**3
wwm_ardhuin_old.F90:1019:      !WRITE(DBG%FHNDL,*) 'TAUPX', TAUPX, TAUPY, TTAUWSHELTER, XSTRESS, YSTRESS
wwm_ardhuin_old.F90:1037:         !WRITE(DBG%FHNDL,*) XK, I, ILEVTAIL, CONST0, TEMP, DELTAIL 
wwm_ardhuin_old.F90:1040:         !WRITE(DBG%FHNDL,*) DELK1, DELK2, XK, I, J
wwm_ardhuin_old.F90:1054:      !WRITE(DBG%FHNDL,*) 'STRESSES', TAUHF, TAUWX, TAUWY
wwm_ardhuin_old.F90:1170:!      WRITE(6,*) 'INSIN4:',FLTABS, SSDSDTH, SSDSC3, SSDSBCK
wwm_ardhuin_old.F90:1208:!      WRITE(6,*) 'INSIN4b:',DIKCUMUL                               
wwm_ardhuin_old.F90:1870:            !WRITE(994,*) I,ABR,fsubw
wwm_ardhuin_old.F90:1975:!      write(DBG%FHNDL,*) z0, ustar, windspeed
wwm_ardhuin_old.F90:2363:!/DEBUG WRITE(6,*) 'ERROR IN SDS',DEPTH,USTAR, IX, IY
wwm_aux.F90:28:                  WRITE(STAT%FHNDL,*) 'Gradients of depth '
wwm_aux.F90:29:                  WRITE(STAT%FHNDL,*) '@D/@X '
wwm_aux.F90:31:                     WRITE(STAT%FHNDL,'(1X,I5,3F10.5)') IP, DDEP(IP,1)
wwm_aux.F90:56:                        WRITE(STAT%FHNDL,*) IP, SLMAX, GDL, GDD , 'MAXSLOPE'
wwm_aux.F90:103:                  WRITE(STAT%FHNDL,*) 'Gradients of depth and current'
wwm_aux.F90:104:                  WRITE(STAT%FHNDL,*) '@U/@X     @V/@X'
wwm_aux.F90:106:                     WRITE(STAT%FHNDL,'(1X,I5,3F10.5)') IP, DCUX, DCUY
wwm_aux.F90:122:                  WRITE(STAT%FHNDL,*) 'The Gradient of Depth and Current'
wwm_aux.F90:123:                  WRITE(STAT%FHNDL,*) ' @U/@X    @U/@Y    @V/@X    @V/@Y'
wwm_aux.F90:125:                     WRITE(STAT%FHNDL,'(1X,I5,4F15.7)') IP, DCUX(IP,1), DCUX(IP,2), DCUY(IP,1), DCUY(IP,2)
wwm_aux.F90:242:           WRITE(2305) 1.
wwm_aux.F90:243:           WRITE(2305) (DVDX(IP), DVDY(IP), SQRT(DVDY(IP)**2+DVDY(IP)**2), IP = 1, MNP)
wwm_aux.F90:532:            !write(dbg%fhndl,*) 'boundary -------1------', TIME, IP, IP_IS_STEADY(IP)
wwm_aux.F90:547:              !write(dbg%fhndl,*) 'converged -------2------', TIME, IP, IP_IS_STEADY(IP)
wwm_aux.F90:550:              !write(dbg%fhndl,*) 'not converged -------3------', TIME, IP, IP_IS_STEADY(IP)
wwm_aux.F90:559:        !write(*,*) time, maxval(IP_IS_STEADY), minval(IP_IS_STEADY)
wwm_aux.F90:568:          !WRITE(*,*) IE, IE_IS_STEADY(IE)
wwm_aux.F90:578:            !IF (IP_IS_STEADY(IP) .GT. 2) WRITE(*,*) TIME, IP, IP_IS_STEADY(IP)
wwm_aux.F90:610:           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 1 REACHED IN', CONV1, '% GRIDPOINTS'
wwm_aux.F90:611:           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 2 REACHED IN', CONV2, '% GRIDPOINTS'
wwm_aux.F90:612:           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 3 REACHED IN', CONV3, '% GRIDPOINTS'
wwm_aux.F90:613:           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 4 REACHED IN', CONV4, '% GRIDPOINTS'
wwm_aux.F90:614:           WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 5 REACHED IN', CONV5, '% GRIDPOINTS'
wwm_aux.F90:617:         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 1 REACHED IN', CONV1, '% GRIDPOINTS'
wwm_aux.F90:618:         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 2 REACHED IN', CONV2, '% GRIDPOINTS'
wwm_aux.F90:619:         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 3 REACHED IN', CONV3, '% GRIDPOINTS'
wwm_aux.F90:620:         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 4 REACHED IN', CONV4, '% GRIDPOINTS'
wwm_aux.F90:621:         WRITE(STAT%FHNDL,*) 'CONVERGENCE CRIT. 5 REACHED IN', CONV5, '% GRIDPOINTS'
wwm_aux.F90:943:        WRITE(4001)  SNGL(RTIME)
wwm_aux.F90:944:        WRITE(4001) (1., 1., SNGL(INIT(IP)), IP = 1, MNP)
wwm_aux.F90:1115:        !write(*,'(3I10,5F15.8)') j, ie, ip, wi(j), DEP(IP), WATLEV(IP), CURTXY(IP,:) 
wwm_aux.F90:1179:             WRITE(STAT%FHNDL,'("+TRACE...",2A)') 'Reading the header of the serial file ... HEADER ', TRIM(HEADLN)
wwm_aux.F90:1186:             WRITE(STAT%FHNDL,'("+TRACE...",2A)') 'Reading the file of the request ... ', TRIM(FILEN)
wwm_aux.F90:1272:           WRITE(2222,'(I10,2F20.8,F15.4)') I-1, XP(I), YP(I), DEP(I)
wwm_aux.F90:1278:           WRITE(2222,'(5I10)') INE(1,I)-1, INE(2,I)-1, INE(3,I)-1, 0, I-1
wwm_aux.F90:1305:           WRITE(2222,'(I10,2F20.8,F15.4)') I-1, XP(I), YP(I), DEP(I)
wwm_aux.F90:1311:           WRITE(2222,'(5I10)') INE(1,I)-1, INE(2,I)-1, INE(3,I)-1, 0, I-1
wwm_aux.F90:1326:      WRITE(IFILE,'(A)') 'C system.dat, made by tri2sys'
wwm_aux.F90:1327:      WRITE(IFILE,'(A)') 'C Number of Boundary Nodes:'
wwm_aux.F90:1328:      WRITE(IFILE,'(I10)') NKR
wwm_aux.F90:1329:      WRITE(IFILE,'(A)') 'C Number of Domain Nodes:'
wwm_aux.F90:1330:      WRITE(IFILE,'(I10)') NKG
wwm_aux.F90:1331:      WRITE(IFILE,'(A)') 'C Koordinaten und Skalarwerte der Knoten'
wwm_aux.F90:1332:      WRITE(IFILE,'(A)') 'C --------------------------------------'
wwm_aux.F90:1333:      WRITE(IFILE,'(A)') 'C Zuerst die Randknoten  (Anzahl s.o.),'
wwm_aux.F90:1334:      WRITE(IFILE,'(A)') 'C dann die Gebietsknoten (Anzahl s.o.).'
wwm_aux.F90:1335:      WRITE(IFILE,'(A)') 'C ------------+-------------+-------------+---------------'
wwm_aux.F90:1336:      WRITE(IFILE,'(A)') 'C     Nr.     |  x-Koord.   |   y-Koord.  | Skalarwert'
wwm_aux.F90:1337:      WRITE(IFILE,'(A)') 'C ------------+-------------+-------------+---------------'
wwm_aux.F90:1348:      WRITE(IFILE,'(A)') "C ------------------------------------------------------------"
wwm_aux.F90:1349:      WRITE(IFILE,'(A)') "C Anzahl der Elemente:"
wwm_aux.F90:1350:      WRITE(IFILE,'(I11)') NELEM
wwm_aux.F90:1351:      WRITE(IFILE,'(A)') "C Elementverzeichnis"
wwm_aux.F90:1352:      WRITE(IFILE,'(A)') "C ------------------------------------------------------------"
wwm_aux.F90:1353:      WRITE(IFILE,'(A)') "C    Knoten i  Knoten j  Knoten k   Kennung     Nr."
wwm_aux.F90:1832:      WRITE(DBG%FHNDL, *) TRIM(string)
wwm_aux.F90:1855:        WRITE(ErrMsg,10) TRIM(string1), TRIM(string2)
wwm_aux.F90:2421:      WRITE(eStrNb,*) TheNb
wwm_aux.F90:2436:      WRITE(eStr,40) TRIM(eStrZero),TRIM(eStrNb)
wwm_aux.F90:2462:!                   WRITE(*,'(3I10,4E15.4)') IP, IS, ID, STAT2D(IS,ID), AC2(IP,IS,ID)/ETOT, AC2(IP,IS,ID), ETOT
wwm_aux.F90:2480:!              WRITE(*,'(2I10,F15.4)') IS, ID, STAT2D(IS,ID)
wwm_aux.F90:2482:!             WRITE(*,*) IS, STAT1D(IS)
wwm_aux_parall.F90:14:      WRITE(STAT%FHNDL,*) 'Before the loop of send/recv stat=', iresult
wwm_aux_parall.F90:28:      WRITE(STAT%FHNDL,*) 'After the loop of send/recv stat=', iresult
wwm_babanin.F90:66:       write(DBG%FHNDL,*)'a1,a2 = ',a1,a2
wwm_babanin.F90:78:!      WRITE(DBG%FHNDL,*) IS, EDENS(IS), EDENST(IS), DEDENS(IS)
wwm_babanin.F90:140:          write(DBG%FHNDL,205)f(is),EDENS(is),EDENST(is),T1(is),T2(is)
wwm_babanin.F90:172:       WRITE(DBG%FHNDL,*)'integral of T1,T2,Sds = ',ST1_INT,ST2_INT,Sds_INT
wwm_babanin.F90:373:       write(DBG%FHNDL,206)tau_total,tauv,tau_normal,TAUWLIM
wwm_babanin.F90:383:          write(DBG%FHNDL,*)'reducing Sin to make tau_normal match TAUWLIM'
wwm_babanin.F90:388:!      write(DBG%FHNDL,*) ip, wind10, ufric, tau_normal 
wwm_babanin.F90:421:          write(DBG%FHNDL,207)iter,REDUC,RCHANGE,tau_normal,err
wwm_babanin.F90:436:       write(DBG%FHNDL,*)'no solution found at gridpoint', IP, SUM(AC2(IP,:,:))
wwm_babanin.F90:437:       write(DBG%FHNDL,'(A20,F15.8,A20,F15.8,A10,F15.8,A20,F15.8)')'tau_normal = ',tau_normal,' TAUWLIM =', TAUWLIM,'err =',err,'(abs(err)/TAUWLIM)  = ',abs(err)/TAUWLIM
wwm_babanin.F90:438:       write(DBG%FHNDL,*) wind10, ufric(ip)
wwm_babanin.F90:441:          write(401)nf_new     ! scalar
wwm_babanin.F90:442:          write(401)Lfactor_L  ! LFACTOR_L(nf_new)
wwm_babanin.F90:443:          write(401)S_in1D_L   ! S_in1D_L(nf_new)
wwm_babanin.F90:444:          write(401)df         ! DF(nf_new)
wwm_babanin.F90:445:          write(401)CINV_L     ! CINV_L(nf_new)
wwm_babanin.F90:446:          write(401)TAUWLIM    ! scalar
wwm_babanin.F90:447:          write(401)tau_normal ! scalar
wwm_babanin.F90:448:          write(401)REDUC      ! scalar
wwm_babanin.F90:449:          write(401)U10MOD     ! scalar
wwm_babanin.F90:450:          write(401)rhow       ! scalar
wwm_babanin.F90:567:!                write(DBG%FHNDL,*)'myu(',IS,') = ',myu
wwm_babanin.F90:569:!                write(DBG%FHNDL,*)'SWDIS(IS),SWDIS_CHECK = ',
wwm_bdcons.F90:97:            WRITE(IOBPOUT%FHNDL,*) IP, ID, IOBWB(IP)
wwm_bdcons.F90:101:              WRITE(IOBPDOUT%FHNDL,*) IP, ID, SPDIR(ID)*RADDEG, IOBPD(ID,IP), IOBP(IP)
wwm_bdcons.F90:131:         WRITE(STAT%FHNDL,*) 'BOUNDARY FILE NAME'
wwm_bdcons.F90:132:         WRITE(STAT%FHNDL,*) 'IGRIDTYPE=', IGRIDTYPE
wwm_bdcons.F90:133:         WRITE(STAT%FHNDL,*) BND%FHNDL, BND%FNAME
wwm_bdcons.F90:156:         IF (myrank == 0) WRITE(STAT%FHNDL,*) 'reading in the boundary flags'
wwm_bdcons.F90:172:         WRITE(STAT%FHNDL,*) 'FINISHED READING BOUNDARY FILE NAME'
wwm_bdcons.F90:180:            WRITE(wwmerr, *) 'IOBP(IP) must not be .gt. 4', IP, ' iobp=', iobp(IP)
wwm_bdcons.F90:203:        WRITE(STAT%FHNDL,*) 'FINISHED WITH EXCHANGE OF BOUNDARY MAPPINGS' 
wwm_bdcons.F90:301:         WRITE(STAT%FHNDL,*)'FINISHED SETTING THE FLAGS'
wwm_bdcons.F90:302:         WRITE(STAT%FHNDL,*)'Gloabl bnd list from init.:', IWBMNPGL,IWBNDGL(:)
wwm_bdcons.F90:303:         WRITE(STAT%FHNDL,*)'Local bnd list from init.:', IWBMNP,IWBNDLC(:)
wwm_bdcons.F90:345:         WRITE(STAT%FHNDL,'("+TRACE...",A,I10)') 'Number of Active Wave Boundary Nodes', IWBMNP
wwm_bdcons.F90:352:            WRITE(IOBPOUT%FHNDL,*) IP, IOBP(IP)
wwm_bdcons.F90:471:      WRITE(STAT%FHNDL,*) 'IGRIDTYPE=', IGRIDTYPE
wwm_bdcons.F90:472:      WRITE(STAT%FHNDL,*) 'BND%FHNDL=', BND%FHNDL
wwm_bdcons.F90:473:      WRITE(STAT%FHNDL,*) 'BND%FNAME=', TRIM(BND%FNAME)
wwm_bdcons.F90:526:          WRITE(wwmerr, *) 'NextGen: We need iobp<=2 but ip=', IP, ' iobp=', iobp(IP)
wwm_bdcons.F90:691:            WRITE(IOBPOUT%FHNDL,*) IP, IOBP(IP)
wwm_bdcons.F90:782:            WRITE(DBG%FHNDL,*) 'BND%FNAME=', BND%FNAME
wwm_bdcons.F90:783:            Write(errmsg,*) 'Error in bnd file', BND%FNAME
wwm_bdcons.F90:861:        WRITE(STAT%FHNDL,'("+TRACE...",A,I10)') 'Number of Active Wave Boundary Nodes', IWBMNP
wwm_bdcons.F90:870:            WRITE(IOBPOUT%FHNDL,*) IP, IOBP(IP)
wwm_bdcons.F90:932:         WRITE(STAT%FHNDL,*) 'WAVE BOUNDARY CONDITION CALLED', IFILE, IT, CALLFROM
wwm_bdcons.F90:942:           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'Parametric Wave Boundary Condition is prescribed'
wwm_bdcons.F90:1021:                 WRITE(STAT%FHNDL,*)'GETWW3SPECTRA CALLED'
wwm_bdcons.F90:1023:                 WRITE(STAT%FHNDL,*)'GETWW3SPECTRA SUCCEEDED'
wwm_bdcons.F90:1025:                   WRITE(DBG%FHNDL,*) ' AFTER CALL GET_BINARY_WW3_SPECTRA',  SUM(WBACOUT)
wwm_bdcons.F90:1034:               WRITE(STAT%FHNDL,'("+TRACE...",A)') '1d Spectra is given as Wave Boundary Condition'
wwm_bdcons.F90:1038:               WRITE(STAT%FHNDL,'("+TRACE...",A)') '2d Spectra is given as Wave Boundary Condition'
wwm_bdcons.F90:1043:                 WRITE(STAT%FHNDL,*)'GETWW3SPECTRA CALLED SUCCEED'
wwm_bdcons.F90:1174:        WRITE(STAT%FHNDL,*) 'IS=', IS, eSum
wwm_bdcons.F90:1181:        WRITE(STAT%FHNDL,*) 'ID=', ID, eSum, DiffAng
wwm_bdcons.F90:1214:        WRITE(DBG%FHNDL,*) 'HS    PER    DIR    DPSR    SHAPE   DEGEXP    GAUSS   PEAK'
wwm_bdcons.F90:1215:        WRITE(DBG%FHNDL,'(8F10.4)') SPPAR(8)
wwm_bdcons.F90:1303:          IF (LDEBUG) WRITE(DBG%FHNDL,*) 'IS LOOP', IS, SF, FPK, SYF, RA
wwm_bdcons.F90:1351:!        write(*,'(I10,4F15.8)') ITPER, MPER, &
wwm_bdcons.F90:1369:!      write(*,*) adir, SPPAR(3), DEG, LNAUTIN
wwm_bdcons.F90:1387:        !write(*,*) aacos, spdir(id), adir
wwm_bdcons.F90:1393:              !WRITE(*,*) 'ERROR', AACOS, COS(DDIR) 
wwm_bdcons.F90:1402:          !write(*,'(2I10,2F15.8)') is, id, cdir, ACLOC(IS,MDC)
wwm_bdcons.F90:1840:           WRITE(DBG%FHNDL,*) ' ENTERING SET BOUNDARY CONDITION ',  SUM(AC2)
wwm_bdcons.F90:1849:                 !WRITE(*,*) DTMP, MAIN%BMJD, BND_TIME_ALL_FILES(1,1)
wwm_bdcons.F90:1865:                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 1', IFILE, IT, LBINTER
wwm_bdcons.F90:1869:                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 2', IFILE, IT, LBINTER
wwm_bdcons.F90:1878:                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 3', IFILE, IT, LBINTER
wwm_bdcons.F90:1882:                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 4', IFILE, IT, LBINTER
wwm_bdcons.F90:1898:                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 1', 1, IT, LBINTER
wwm_bdcons.F90:1902:                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 2', 1, IT, LBINTER
wwm_bdcons.F90:1910:                   WRITE(DBG%FHNDL,*) ' AFTER CALL TO WAVE_BOUNDARY_CONDITION LBINTER TRUE',  SUM(WBAC), SUM(WBACOLD), SUM(WBACNEW)
wwm_bdcons.F90:1915:                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 3', 1, IT, LBINTER
wwm_bdcons.F90:1919:                   WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION 4', 1, IT, LBINTER
wwm_bdcons.F90:1924:                   WRITE(DBG%FHNDL,*) ' AFTER CALL TO WAVE_BOUNDARY_CONDITION LBINTER FALSE',  SUM(WBAC), SUM(WBACOLD), SUM(WBACNEW)
wwm_bdcons.F90:1939:                 WRITE(DBG%FHNDL,*) ' AFTER TIME INTERPOLATION NO READ OF FILE',  SUM(WBAC), SUM(DSPEC)
wwm_bdcons.F90:1958:             WRITE(DBG%FHNDL,*) ' FINISHED WITH BOUNDARY CONDITION ',  SUM(AC2)
wwm_breaking.F90:98:        !IF (QB .GT. 0.1) WRITE(*,'(7F15.4)') HMAX(IP), BETA, QB, KME, HS
wwm_breaking.F90:121:              !write(*,'(5F15.10)') ALPBJ * QB * SME / BETA2, ALPBJ, QB, SME, BETA2
wwm_breaking.F90:136:              !if (abs(surfa0) .gt. zero) write(*,*) surfa0, surfa1
wwm_breaking.F90:141:              !if (surfa0 .lt. zero) write(*,*) is, id, SURFA0
wwm_breaking.F90:242:        !IF (QB .GT. 0.1) WRITE(*,'(7F15.4)') HMAX(IP), BETA, QB, KME, HS
wwm_breaking.F90:264:        !IF (QB .GT. 0.001) WRITE(*,*) QB, SURFA0, SURFA1 
wwm_buildstress.F90:89:      IF (ITEST.GE.1) WRITE(IU06,*) ' SUB. BUILDSTRESS: INPUT OF RESTART FILES DONE'
wwm_compute.F90:17:         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START COMPUTE COMPUTE_SIMPLE_EXPLICIT'
wwm_compute.F90:22:           WRITE(DBG%FHNDL,*) ' AFTER ENTERING COMPUTE ',  SUM(AC2)
wwm_compute.F90:52:           WRITE(DBG%FHNDL,*) ' AFTER DIRECTION AND FREQUENCY -1- ',  SUM(AC2)
wwm_compute.F90:63:           WRITE(DBG%FHNDL,*) ' AFTER DIRECTION AND FREQUENCY -2- ',  SUM(AC2)
wwm_compute.F90:73:           WRITE(DBG%FHNDL,*) ' AFTER SPATIAL ',  SUM(AC2)
wwm_compute.F90:97:               WRITE(111112,'(A10,I10)') 'AFTER', IP
wwm_compute.F90:98:               WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(1,:,:))
wwm_compute.F90:99:               WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(1,:,:))
wwm_compute.F90:100:               WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
wwm_compute.F90:101:               WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
wwm_compute.F90:102:               WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
wwm_compute.F90:103:               WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
wwm_compute.F90:104:               WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
wwm_compute.F90:105:               WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
wwm_compute.F90:106:               WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
wwm_compute.F90:107:               WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
wwm_compute.F90:108:               WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
wwm_compute.F90:109:               WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
wwm_compute.F90:110:               WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
wwm_compute.F90:111:               WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(1,:,:))
wwm_compute.F90:112:               WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(1,:))
wwm_compute.F90:131:               WRITE(111112,'(A10,I10)') 'AFTER', IP
wwm_compute.F90:132:               WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(1,:,:))
wwm_compute.F90:133:               WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(1,:,:))
wwm_compute.F90:134:               WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
wwm_compute.F90:135:               WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
wwm_compute.F90:136:               WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
wwm_compute.F90:137:               WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
wwm_compute.F90:138:               WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
wwm_compute.F90:139:               WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
wwm_compute.F90:140:               WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
wwm_compute.F90:141:               WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
wwm_compute.F90:142:               WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
wwm_compute.F90:143:               WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
wwm_compute.F90:144:               WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
wwm_compute.F90:145:               WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(1,:,:))
wwm_compute.F90:146:               WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(1,:))
wwm_compute.F90:173:           WRITE(DBG%FHNDL,*) ' AFTER SOURCES ',  SUM(AC2)
wwm_compute.F90:181:           WRITE(DBG%FHNDL,*) ' AFTER BREAK LIMIT ',  SUM(AC2)
wwm_compute.F90:185:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----SIMPLE SPLITTING SCHEME-----'
wwm_compute.F90:186:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1
wwm_compute.F90:187:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SOURCES                          ', TIME6-TIME5
wwm_compute.F90:188:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME5-TIME4
wwm_compute.F90:189:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS THETA SPACE          ', TIME4-TIME3
wwm_compute.F90:190:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SIGMA SPACE          ', TIME3-TIME2
wwm_compute.F90:191:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU MICHE LIMITER                ', TIME6-TIME5
wwm_compute.F90:192:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME6-TIME1
wwm_compute.F90:193:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
wwm_compute.F90:195:         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE COMPUTE_SIMPLE_EXPLICIT'
wwm_compute.F90:213:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE COMPUTE_SEMI_IMPLICIT'
wwm_compute.F90:271:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----IMPLICIT SPLITTING SCHEME-----'
wwm_compute.F90:272:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1
wwm_compute.F90:273:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME6-TIME5
wwm_compute.F90:274:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SPECTRAL SPACE       ', TIME3-TIME2
wwm_compute.F90:275:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SOURCES              ', TIME5-TIME4
wwm_compute.F90:276:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'ACTION LIMITER                   ', TIME7-TIME6
wwm_compute.F90:277:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'MICHE LIMITER                    ', TIME8-TIME7+TIME4-TIME3
wwm_compute.F90:278:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME8-TIME1
wwm_compute.F90:279:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
wwm_compute.F90:280:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE COMPUTE_SEMI_IMPLICIT'
wwm_compute.F90:291:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE COMPUTE_SPATIAL'
wwm_compute.F90:307:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SPATIAL'
wwm_compute.F90:318:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_SOURCES_EXP'
wwm_compute.F90:325:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SOURCES_EXP'
wwm_compute.F90:446:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----IMPLICIT -----'
wwm_compute.F90:447:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1
wwm_compute.F90:448:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS IMPLICIT             ', TIME5-TIME4
wwm_compute.F90:449:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SOURCES              ', TIME4-TIME3+TIME6-TIME5
wwm_compute.F90:450:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'ACTION LIMITER                   ', TIME7-TIME6
wwm_compute.F90:451:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'MICHE LIMITER                    ', TIME3-TIME2
wwm_compute.F90:452:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME7-TIME1
wwm_compute.F90:453:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
wwm_compute.F90:454:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE COMPUTE_SEMI_IMPLICIT'
wwm_coupl_roms.F90:16:      WRITE(DBG%FHNDL,'("+TRACE...",A)') 'OPEN PIPE ROMS'
wwm_coupl_roms.F90:20:      WRITE(DBG%FHNDL,*) 'WWM: open pipe ExchImport'
wwm_coupl_roms.F90:24:      WRITE(DBG%FHNDL,*) 'WWM: open pipe ExchExport'
wwm_coupl_roms.F90:26:      WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END OPEN PIPE ROMS'
wwm_coupl_roms.F90:57:        WRITE(DBG%FHNDL,'("+TRACE...",A)') 'READING PIPE'
wwm_coupl_roms.F90:111:        WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END READING PIPE'
wwm_coupl_roms.F90:159:          WRITE(101)  HS, HSWE,                                        &
wwm_coupl_roms.F90:202:            WRITE(101) OUTT_TOT(IP, 1), OUTT_TOT(IP, 2),                 &
wwm_coupl_roms.F90:216:      WRITE(DBG%FHNDL,*) 'export WWM: ending of writing data'
wwm_coupl_roms.F90:346:        WRITE(DBG%FHNDL,*) 'MinValIndex(Dir,Inv)=', MinValIndex, MinValIndexInv
wwm_coupl_roms.F90:398:          WRITE(DBG%FHNDL,*) 'nb1=', nb1, ' ne_global=', ne_global
wwm_coupl_roms.F90:399:          WRITE(DBG%FHNDL,*) 'nb2=', nb2, ' np_global=', np_global
wwm_coupl_roms.F90:447:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1, rnk=', myrank
wwm_coupl_roms.F90:452:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1.2, rnk=', myrank
wwm_coupl_roms.F90:457:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1.3, rnk=', myrank
wwm_coupl_roms.F90:462:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 2, rnk=', myrank
wwm_coupl_roms.F90:473:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 3, rnk=', myrank
wwm_coupl_roms.F90:479:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 4, rnk=', myrank
wwm_coupl_roms.F90:485:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 5, rnk=', myrank
wwm_coupl_roms.F90:491:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 6, rnk=', myrank
wwm_coupl_roms.F90:497:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 7, rnk=', myrank
wwm_coupl_roms.F90:503:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 8, rnk=', myrank
wwm_coupl_roms.F90:509:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 9, rnk=', myrank
wwm_coupl_roms.F90:518:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 10, rnk=', myrank
wwm_coupl_roms.F90:526:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 11, rnk=', myrank
wwm_coupl_roms.F90:537:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 12, rnk=', myrank
wwm_coupl_roms.F90:548:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 13, rnk=', myrank
wwm_coupl_roms.F90:557:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 14', myrank
wwm_coupl_roms.F90:570:        WRITE(DBG%FHNDL,*) 'nbNeedTot=', TheArr_OCNtoWAV_rho % nbNeedTot
wwm_coupl_roms.F90:571:        WRITE(DBG%FHNDL,*) 'nbProc=', TheArr_OCNtoWAV_rho % nbProc
wwm_coupl_roms.F90:572:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, WAV, step 14'
wwm_coupl_roms.F90:588:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 15, rnk=', myrank
wwm_coupl_roms.F90:604:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 16, rnk=', myrank
wwm_coupl_roms.F90:612:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 17, rnk=', myrank
wwm_coupl_roms.F90:621:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 18, rnk=', myrank
wwm_coupl_roms.F90:632:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 19, rnk=', myrank
wwm_coupl_roms.F90:641:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 20, rnk=', myrank
wwm_coupl_roms.F90:652:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 21, rnk=', myrank
wwm_coupl_roms.F90:668:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 22, rnk=', myrank
wwm_coupl_roms.F90:684:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 23, rnk=', myrank
wwm_coupl_roms.F90:689:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 23.1'
wwm_coupl_roms.F90:694:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 23.2'
wwm_coupl_roms.F90:699:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 23.3'
wwm_coupl_roms.F90:717:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 24, rnk=', myrank
wwm_coupl_roms.F90:727:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 25, rnk=', myrank
wwm_coupl_roms.F90:728:        WRITE(DBG%FHNDL,*) 'MyRankGlobal=', MyRankGlobal
wwm_coupl_roms.F90:743:        WRITE(DBG%FHNDL,*) 'SumDepReceive=', SumDepReceive
wwm_coupl_roms.F90:744:        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, WAV, step 33'
wwm_coupl_roms.F90:745:        WRITE(DBG%FHNDL,*) 'WAV, rnk=', myrank
wwm_coupl_roms.F90:761:        WRITE(DBG%FHNDL,*) 'dep_rho, min=', minBathy, ' max=', maxBathy
wwm_coupl_roms.F90:773:        WRITE(DBG%FHNDL,*) 'DEP, min=', minBathy, ' max=', maxBathy
wwm_coupl_roms.F90:794:            WRITE(DBG%FHNDL,*) 'AD, IP=', IP, dep_rho(IP), DEP(IP)
wwm_coupl_roms.F90:795:            WRITE(DBG%FHNDL,*) 'AD, xp, yp=', XP(IP), YP(IP)
wwm_coupl_roms.F90:801:!        WRITE(DBG%FHNDL,*) 'AD, AbsDiff=', AbsDiff
wwm_coupl_roms.F90:802:!        WRITE(DBG%FHNDL,*) 'AD, IP=', iNodeSel, dep_rho(iNodeSel), DEP(iNodeSel)
wwm_coupl_roms.F90:803:!        WRITE(DBG%FHNDL,*) 'AD, xp, yp=', XP(iNodeSel), YP(iNodeSel)
wwm_coupl_roms.F90:804:        WRITE(DBG%FHNDL,*) 'AD, SumDep1=', SumDep1, ' SumDep2=', SumDep2
wwm_coupl_roms.F90:805:        WRITE(DBG%FHNDL,*) 'AD, SumDiff=', SumDiff
wwm_coupl_roms.F90:821:        WRITE(DBG%FHNDL,*) 'End ROMS_COUPL_INITIALIZE'
wwm_coupl_roms.F90:1084:        WRITE(DBG%FHNDL,*) 'WWM: Begin PGMCL_ROMS_IN'
wwm_coupl_roms.F90:1101:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, After Data receive'
wwm_coupl_roms.F90:1139:        WRITE(DBG%FHNDL,*) 'WAV, MaxUwind=', MaxUwind, ' avgUwind=', avgUwind
wwm_coupl_roms.F90:1140:        WRITE(DBG%FHNDL,*) 'WAV, MaxVwind=', MaxVwind, ' avgVwind=', avgVwind
wwm_coupl_roms.F90:1141:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 2'
wwm_coupl_roms.F90:1158:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 3'
wwm_coupl_roms.F90:1168:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 4'
wwm_coupl_roms.F90:1193:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 5'
wwm_coupl_roms.F90:1214:        WRITE(DBG%FHNDL,*) 'After the receive'
wwm_coupl_roms.F90:1222:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 6'
wwm_coupl_roms.F90:1242:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 7'
wwm_coupl_roms.F90:1282:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 1'
wwm_coupl_roms.F90:1321:        WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'MaxNormTau=', MaxNormTau
wwm_coupl_roms.F90:1322:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 5.1'
wwm_coupl_roms.F90:1361:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 5.3'
wwm_coupl_roms.F90:1415:            WRITE(DBG%FHNDL,*) 'eStokesNorm=', eStokesNorm
wwm_coupl_roms.F90:1416:            WRITE(DBG%FHNDL,*) 'KLM=', KLM, 'WLM=', WLM
wwm_coupl_roms.F90:1417:            WRITE(DBG%FHNDL,*) 'cPhase=', cPhase, 'kD=', kD
wwm_coupl_roms.F90:1418:            WRITE(DBG%FHNDL,*) 'HS=', HS, ' DEP=', DEP(IP)
wwm_coupl_roms.F90:1452:        WRITE(DBG%FHNDL,*) 'WAV, MaxHwave=', MaxHwave, ' avgHwave=', avgHwave
wwm_coupl_roms.F90:1453:        WRITE(DBG%FHNDL,*) 'WAV, MaxLwave=', MaxLwave, ' avgLwave=', avgLwave
wwm_coupl_roms.F90:1454:        WRITE(DBG%FHNDL,*) 'WAV, MaxStokesNorm=', MaxStokesNorm
wwm_coupl_roms.F90:1455:        WRITE(DBG%FHNDL,*) 'WAV, avgStokesNorm=', avgStokesNorm
wwm_coupl_roms.F90:1456:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 6'
wwm_coupl_roms.F90:1463:        WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'AvgNormFV2=', AvgUFRICsqr
wwm_coupl_roms.F90:1464:        WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'AvgCdU2=', AvgStressCd
wwm_coupl_roms.F90:1465:        WRITE(DBG%FHNDL,*) 'AvgCd=', AvgCd, ' AvgAlpha=', AvgAlpha
wwm_coupl_roms.F90:1466:        WRITE(DBG%FHNDL,*) 'AvgWind=', AvgWind
wwm_coupl_roms.F90:1483:        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 11'
wwm_coupl_selfe.F90:239:                  write(errmsg,*)'R.S.: div by 0 (0);',iplg(IP),EWK(IP),DEP(IP)
wwm_coupl_selfe.F90:403:!         write(12,*)'Checking Sxx,Sxy,Syy:'
wwm_coupl_selfe.F90:451:!          write(12,*)'Checking R.S.'
wwm_coupl_selfe.F90:549:      !write(*,*) sum(forcexy)
wwm_coupl_shyfem.F90:16:           WRITE(DBG%FHNDL,'("+TRACE...",A)') 'OPEN PIPE'
wwm_coupl_shyfem.F90:51:           WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END OPEN PIPE'
wwm_coupl_shyfem.F90:99:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'READING PIPE'
wwm_coupl_shyfem.F90:163:        WRITE(STAT%FHNDL,*) 'CHECK MAX UX,UY,H'
wwm_coupl_shyfem.F90:164:        WRITE(STAT%FHNDL,*) MAXVAL(CURTXY(:,1)), MAXVAL(CURTXY(:,2)), MAXVAL(WATLEV)
wwm_coupl_shyfem.F90:165:        WRITE(STAT%FHNDL,*) 'CHECK MIN UX,UY,H'
wwm_coupl_shyfem.F90:166:        WRITE(STAT%FHNDL,*) MINVAL(CURTXY(:,1)), MINVAL(CURTXY(:,2)), MINVAL(WATLEV)
wwm_coupl_shyfem.F90:167:        WRITE(2001,*) 'CHECK MAX UX,UY,H'
wwm_coupl_shyfem.F90:168:        WRITE(2001,*) MAXVAL(CURTXY(:,1)), MAXVAL(CURTXY(:,2)), MAXVAL(WATLEV)
wwm_coupl_shyfem.F90:169:        WRITE(2001,*) 'CHECK MIN UX,UY,H'
wwm_coupl_shyfem.F90:170:        WRITE(2001,*) MINVAL(CURTXY(:,1)), MINVAL(CURTXY(:,2)), MINVAL(WATLEV)
wwm_coupl_shyfem.F90:171:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'END READ PIPE WWM'
wwm_coupl_shyfem.F90:265:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'WRITING PIPE'
wwm_coupl_shyfem.F90:280:            WRITE(11101)  SXX3D(IL,IP)             !ccf
wwm_coupl_shyfem.F90:281:            WRITE(11102)  SYY3D(IL,IP)             !ccf
wwm_coupl_shyfem.F90:282:            WRITE(11142)  SXY3D(IL,IP)             !ccf
wwm_coupl_shyfem.F90:292:          WRITE(11103) HS
wwm_coupl_shyfem.F90:294:          WRITE(11104) SME01
wwm_coupl_shyfem.F90:296:          WRITE(11105) DM
wwm_coupl_shyfem.F90:298:          !WRITE(11106) KMWAM
wwm_coupl_shyfem.F90:299:          WRITE(11106) TAUW(IP)
wwm_coupl_shyfem.F90:301:          WRITE(11107) TPP
wwm_coupl_shyfem.F90:303:          WRITE(11108) KPP
wwm_coupl_shyfem.F90:305:          WRITE(11109) ORBITAL
wwm_coupl_shyfem.F90:307:          WRITE(11110) STOKES_X(:,IP)
wwm_coupl_shyfem.F90:309:          WRITE(11111) STOKES_Y(:,IP)
wwm_coupl_shyfem.F90:312:            WRITE(11112) WINDXY(IP,1)
wwm_coupl_shyfem.F90:314:            WRITE(11113) WINDXY(IP,2)
wwm_coupl_shyfem.F90:317:          WRITE(11114) CD(IP) 
wwm_coupl_shyfem.F90:319:          WRITE(11115) JPRESS(IP)
wwm_coupl_shyfem.F90:363:              WRITE(11101)  OUTT(IP, IL       )
wwm_coupl_shyfem.F90:364:              WRITE(11102)  OUTT(IP, IL+  NLVT)
wwm_coupl_shyfem.F90:365:              WRITE(11142)  OUTT(IP, IL+2*NLVT)
wwm_coupl_shyfem.F90:370:            WRITE(11103) OUTT(IP, 1 + 5*NLVT)
wwm_coupl_shyfem.F90:372:            WRITE(11104) OUTT(IP, 2 + 5*NLVT)
wwm_coupl_shyfem.F90:374:            WRITE(11105) OUTT(IP, 3 + 5*NLVT)
wwm_coupl_shyfem.F90:376:            WRITE(11106) OUTT(IP, 4 + 5*NLVT)
wwm_coupl_shyfem.F90:378:            WRITE(11107) OUTT(IP, 5 + 5*NLVT)
wwm_coupl_shyfem.F90:380:            WRITE(11108) OUTT(IP, 6 + 5*NLVT)
wwm_coupl_shyfem.F90:382:            WRITE(11109) OUTT(IP, 7 + 5*NLVT)
wwm_coupl_shyfem.F90:388:            WRITE(11110) STOKES_X_ret
wwm_coupl_shyfem.F90:390:            WRITE(11111) STOKES_Y_ret
wwm_coupl_shyfem.F90:393:              WRITE(11112) OUTT(IP, 10 + 5*NLVT)
wwm_coupl_shyfem.F90:395:              WRITE(11113) OUTT(IP, 11 + 5*NLVT)
wwm_coupl_shyfem.F90:398:            WRITE(11114) OUTT(IP, 8 + 5*NLVT)
wwm_coupl_shyfem.F90:400:            WRITE(11115) OUTT(IP, 9 + 5*NLVT)
wwm_coupl_shyfem.F90:407:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'END WRITING PIPE'
wwm_coupl_timor.F90:14:           WRITE(DBG%FHNDL,'("+TRACE...",A)') 'OPEN PIPE'
wwm_coupl_timor.F90:40:           WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END OPEN PIPE'
wwm_coupl_timor.F90:44:           WRITE(1004) MAIN%DTCOUP
wwm_coupl_timor.F90:52:             write(DBG%FHNDL,*) 'TIME STEP OF THE hydraulic flow MODEL CANNOT BE DIVIDIED WITHOUT A REST'
wwm_coupl_timor.F90:53:             write(DBG%FHNDL,*)'dt Stroemung (s) =',MAIN%DTCUR, ',  dt Kopplung (s) = ',MAIN%DTCOUP
wwm_coupl_timor.F90:57:           WRITE(DBG%FHNDL,'("+TRACE... DTCUR and DTCOUP",A)') MAIN%DTCUR, MAIN%DTCOUP
wwm_coupl_timor.F90:100:        WRITE(DBG%FHNDL,'("+TRACE...",A)') 'READING PIPE'
wwm_coupl_timor.F90:101:        WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END READING PIPE'
wwm_cycle3.F90:95:           WRITE(*,'(A20,6E20.10)') 'LINEAR INPUT', SUM(SSINL), MINVAL(SSINL), MAXVAL(SSINL)
wwm_cycle3.F90:96:           WRITE(*,'(A20,6E20.10)') 'WAVE ACTION', SUM(ACLOC), MINVAL(ACLOC), MAXVAL(ACLOC)
wwm_cycle3.F90:97:           WRITE(*,'(A20,6E20.10)') 'EXP INPUT', SUM(SSINE), SUM(DSSINE), MINVAL(SSINE), MAXVAL(SSINE), MINVAL(DSSINE), MAXVAL(DSSINE)
wwm_cycle3.F90:98:           WRITE(*,'(A20,6E20.10)') 'WHITECAP', SUM(SSDS), SUM(DSSDS), MINVAL(SSDS), MAXVAL(SSDS), MINVAL(DSSDS), MAXVAL(DSSDS)
wwm_cycle3.F90:99:           WRITE(*,'(A20,6E20.10)') 'SNL4', SUM(SSNL4), SUM(DSSNL4), MINVAL(SSNL4), MAXVAL(SSNL4), MINVAL(DSSNL4), MAXVAL(DSSNL4)
wwm_cycle3.F90:100:           WRITE(*,'(A20,6E20.10)') 'SNL3', SUM(SSNL3), SUM(DSSNL3), MINVAL(SSNL3), MAXVAL(SSNL3), MINVAL(DSSNL3), MAXVAL(DSSNL3)
wwm_cycle3.F90:101:           WRITE(*,'(A20,6E20.10)') 'BOTTOM FRICTION', SUM(SSBF), SUM(DSSBF), MINVAL(SSBF), MAXVAL(SSBF), MINVAL(DSSBF), MAXVAL(DSSBF)
wwm_cycle3.F90:102:           WRITE(*,'(A20,6E20.10)') 'BREAKING', SUM(SSBR), SUM(DSSBR), MINVAL(SSBR), MAXVAL(SSBR), MINVAL(DSSBR), MAXVAL(DSSBR)
wwm_cycle3.F90:103:           WRITE(*,'(A20,6E20.10)') 'BREAKING LIMITER', SUM(SSBRL), SUM(DSSBRL), MINVAL(SSBRL), MAXVAL(SSBRL), MINVAL(DSSBRL), MAXVAL(DSSBRL)
wwm_cycle3.F90:104:           WRITE(*,'(A20,6E20.10)') 'LIMITER',  SUM(SSLIM), SUM(DSSLIM), MINVAL(SSLIM), MAXVAL(SSLIM), MINVAL(DSSLIM), MAXVAL(DSSLIM)
wwm_cycle3.F90:105:           WRITE(*,'(A20,6E20.10)') 'TOTAL SOURCE TERMS', SUM(IMATRA), SUM(IMATDA), MINVAL(IMATRA), MAXVAL(IMATRA), MINVAL(IMATDA), MAXVAL(IMATDA)
wwm_diclass.F90:243:      WRITE(*,*) ZLEV 
wwm_diclass.F90:246:      WRITE(*,*) MAXVAL(ZLEV), MAXVAL(XRAY), MAXVAL(YRAY)
wwm_diclass.F90:247:      WRITE(*,*) MINVAL(ZLEV), MINVAL(XRAY), MINVAL(YRAY)
wwm_diffrac.F90:103:           WRITE(555) SNGL(RTIME)
wwm_diffrac.F90:104:           WRITE(555) (SNGL(DIFRX(IP)), SNGL(DIFRY(IP)),SNGL(DIFRM(IP))-1., IP = 1, MNP)
wwm_diffrac.F90:107:         !WRITE(WWMDBG%FHNDL,*) MAXVAL(DIFRM), MAXVAL(DIFRX), MAXVAL(DIFRY)
wwm_diffrac.F90:108:         !WRITE(WWMDBG%FHNDL,*) MINVAL(DIFRM), MINVAL(DIFRX), MINVAL(DIFRY)
wwm_diffrac.F90:180:            WRITE(wwmerr,*)'DFBOT is NaN', IP,KH, CURH(IP), BOTFS(KH), BOTFC(KH), EWK(IP), SLPH(IP)
wwm_diffrac.F90:208:          WRITE(wwmerr,*)'BOTFC is NaN Aron', KH, AUX, AUX1, MyTANH(KH)
wwm_diffrac.F90:245:          WRITE(wwmerr,*)'BOTFS is NaN', BOTFS, AUX, AUX1
wwm_diffrac.F90:298:           WRITE(*,*) 'BOTFC2'
wwm_diffrac.F90:299:           WRITE(*,*) SINHKH, COSHKH, SINH2KH, SINH3KH, KH
wwm_diffrac.F90:300:           WRITE(*,*) AUX, AUX1, KH*MyTANH(KH), (TWO*(COSHKH)**2)
wwm_diffrac.F90:338:          WRITE(*,*) 'BOTFS2'
wwm_diffrac.F90:339:          WRITE(*,*) COSHKH, COSH2KH, SINHKH, SINH2KH, SINH3KH
wwm_diffrac.F90:340:          WRITE(*,*) AUX, AUX1, KH, SECH, SECH**2, COSHKH
wwm_diffrac.F90:394:            WRITE(*,*) 'DFBOT'
wwm_diffrac.F90:395:            WRITE(*,*) DFBOT(IP), BOTFC2(KH), BOTFS2(KH), CURH(IP), EWK(IP), SLPH(IP)
wwm_femeanws.F90:106:          !write(*,'(10F20.10)') EM(IJ) , FM(IJ), DFIMOFR(M), DFIM(M), TEMP2(IJ)
wwm_fkmean.F90:102:!      WRITE(111118,'(I10,F30.20)') ISHALLO, SUM(F)
wwm_fkmean.F90:103:!      WRITE(111118,'(5F20.9)')DELT25,COEFM1,COEF1,COEFA,COEFX
wwm_fkmean.F90:127:            !WRITE(111118,'(4F20.10)') EM(IJ), FM1(IJ), TEMP2(IJ)
wwm_fkmean.F90:158:            IF (LOUTWAM) WRITE(111118,'(4F20.10)') DFIM(M), DFIMOFR(M), DFFR(M)
wwm_fkmean.F90:193:        IF (LOUTWAM)  WRITE(111118,'(4F20.10)') XK(IJ), AK(IJ), F1(IJ), EM(IJ)
wwm_fkmean.F90:199:      IF (LOUTWAM)  WRITE(111118,'(4F20.10)') XK(IJ), AK(IJ), F1(IJ), EM(IJ)
wwm_fluctsplit.F90:295:!              WRITE(DBG%FHNDL,*) '1st IE LOOP CYCLE', IE, IE_IS_STEADY(IE)
wwm_fluctsplit.F90:360:             !WRITE(DBG%FHNDL,'(I10,3F15.6)') IP, SI(IP), KKSUM(IP), DEPTH(IP), DTMAX_EXP
wwm_fluctsplit.F90:369:           !WRITE(STAT%FHNDL,'(2I10,2F15.4)') IS, ID, DTMAX_GLOBAL_EXP, DT4A/DTMAX_GLOBAL_EXP
wwm_fluctsplit.F90:385:             !WRITE(22227,*) IP, CCON(IP), SI(IP)
wwm_fluctsplit.F90:386:             !IF (IP == 24227 .AND. IS == 1) WRITE(DBG%FHNDL,'(2I10,6F20.8)') IP, ID, XP(IP), YP(IP), SI(IP), KKSUM(IP), DEP(IP), CFLCXY(3,IP) 
wwm_fluctsplit.F90:388:           !WRITE(STAT%FHNDL,'(2I10,2F15.4)') IS, ID, DTMAX_GLOBAL_EXP, DT4A/DTMAX_GLOBAL_EXP !AR: Makes very strange error in the code ...
wwm_fluctsplit.F90:413:!         WRITE(STAT%FHNDL,'(3I10,4F15.4)') IS, ID, ITER_EXP(IS,ID), SQRT(MAXVAL(C(1,:))**2+MAXVAL(C(2,:))**2), &
wwm_fluctsplit.F90:421:!                WRITE(DBG%FHNDL,*) '2nd IE LOOP CYCLE', IT, IE, IE_IS_STEADY(IE)
wwm_fluctsplit.F90:431:!                 WRITE(DBG%FHNDL,*) '1st IP LOOP CYCLE', IT, IP, IP_IS_STEADY(IP)
wwm_fluctsplit.F90:436:!             WRITE(*,'(2I10,F20.10,2I20,F20.10)') ID, IS, U(IP_TEST), IOBPD(ID,IP_TEST), IOBDP(IP_TEST), DEP(IP_TEST)
wwm_fluctsplit.F90:479:           WRITE(4001)  SNGL(RTIME)
wwm_fluctsplit.F90:480:           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(IP,1,1)),      &
wwm_fluctsplit.F90:677:           WRITE(4001)  SNGL(RTIME)
wwm_fluctsplit.F90:678:           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(IP,1,1)),      &
wwm_fluctsplit.F90:941:           WRITE(4001)  SNGL(RTIME)
wwm_fluctsplit.F90:942:           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(IP,1,1)),      &
wwm_fluctsplit.F90:1131:!         WRITE(DBG%FHNDL,*) 'CALL SOLVER'
wwm_fluctsplit.F90:1133:!         WRITE(DBG%FHNDL,*) DT4A, MSC, MDC, MNE
wwm_fluctsplit.F90:1134:!         WRITE(DBG%FHNDL,*) 'WRITE CG', SUM(CG)
wwm_fluctsplit.F90:1135:!         WRITE(DBG%FHNDL,*) SUM(XP), SUM(YP)
wwm_fluctsplit.F90:1136:!         WRITE(DBG%FHNDL,*) SUM(IMATRAA), SUM(IMATDAA)
wwm_fluctsplit.F90:1137:!         WRITE(DBG%FHNDL,*) SUM(B), SUM(X)
wwm_fluctsplit.F90:1138:!         WRITE(DBG%FHNDL,*) SUM(IPAR), SUM(FPAR)
wwm_fluctsplit.F90:1139:!         WRITE(DBG%FHNDL,*) SUM(WKSP), SUM(INIU)
wwm_fluctsplit.F90:1140:!         WRITE(DBG%FHNDL,*) SUM(ASPAR), SUM(IA), SUM(JA)
wwm_fluctsplit.F90:1141:!         WRITE(DBG%FHNDL,*) SUM(AU), SUM(FLJAU), SUM(FLJU)
wwm_fluctsplit.F90:1155:!         WRITE(DBG%FHNDL,*) 'SOLUTION'
wwm_fluctsplit.F90:1156:!         WRITE(DBG%FHNDL,*) SUM(B), SUM(X)
wwm_fluctsplit.F90:1157:!         WRITE(DBG%FHNDL,*) SUM(IPAR), SUM(FPAR)
wwm_fluctsplit.F90:1158:!         WRITE(DBG%FHNDL,*) SUM(WKSP), SUM(INIU)
wwm_fluctsplit.F90:1159:!         WRITE(DBG%FHNDL,*) SUM(ASPAR), SUM(JA), SUM(JA)
wwm_fluctsplit.F90:1160:!         WRITE(DBG%FHNDL,*) SUM(AU), SUM(FLJAU), SUM(FLJU)
wwm_fluctsplit.F90:1163:           WRITE(4001)  SNGL(RTIME)
wwm_fluctsplit.F90:1164:           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(IP,1,1)),      &
wwm_fluctsplit.F90:1276:         WRITE(3000+myrank,*) 'IS, ID, sum=', IS, ID, sum(ASPAR)
wwm_fluctsplit.F90:1381:!         WRITE(DBG%FHNDL,*) 'CALL SOLVER'
wwm_fluctsplit.F90:1383:!         WRITE(DBG%FHNDL,*) DT4A, MSC, MDC, MNE
wwm_fluctsplit.F90:1384:!         WRITE(DBG%FHNDL,*) 'WRITE CG', SUM(CG)
wwm_fluctsplit.F90:1385:!         WRITE(DBG%FHNDL,*) SUM(XP), SUM(YP)
wwm_fluctsplit.F90:1386:!         WRITE(DBG%FHNDL,*) SUM(IMATRAA), SUM(IMATDAA)
wwm_fluctsplit.F90:1387:!         WRITE(DBG%FHNDL,*) SUM(B), SUM(X)
wwm_fluctsplit.F90:1388:!         WRITE(DBG%FHNDL,*) SUM(IPAR), SUM(FPAR)
wwm_fluctsplit.F90:1389:!         WRITE(DBG%FHNDL,*) SUM(WKSP), SUM(INIU)
wwm_fluctsplit.F90:1390:!         WRITE(DBG%FHNDL,*) SUM(ASPAR), SUM(IA), SUM(JA)
wwm_fluctsplit.F90:1391:!         WRITE(DBG%FHNDL,*) SUM(AU), SUM(FLJAU), SUM(FLJU)
wwm_fluctsplit.F90:1404:!         WRITE(DBG%FHNDL,*) 'SOLUTION'
wwm_fluctsplit.F90:1405:!         WRITE(DBG%FHNDL,*) SUM(B), SUM(X)
wwm_fluctsplit.F90:1406:!         WRITE(DBG%FHNDL,*) SUM(IPAR), SUM(FPAR)
wwm_fluctsplit.F90:1407:!         WRITE(DBG%FHNDL,*) SUM(WKSP), SUM(INIU)
wwm_fluctsplit.F90:1408:!         WRITE(DBG%FHNDL,*) SUM(ASPAR), SUM(JA), SUM(JA)
wwm_fluctsplit.F90:1409:!         WRITE(DBG%FHNDL,*) SUM(AU), SUM(FLJAU), SUM(FLJU)
wwm_fluctsplit.F90:1412:           WRITE(4001)  SNGL(RTIME)
wwm_fluctsplit.F90:1413:!           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(IP,1,1)),      &
wwm_fluctsplit.F90:1679:           WRITE(4001)  SNGL(RTIME)
wwm_fluctsplit.F90:1680:!           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), SNGL(AC2(IP,1,1)),&
wwm_fluctsplit.F90:2064:         WRITE(STAT%FHNDL,'("+TRACE......",A)') 'CALCULATE CONNECTED AREA SI '
wwm_fluctsplit.F90:2080:         WRITE(STAT%FHNDL,'("+TRACE......",A)') 'MEDIAN DUAL AREA and CCON' 
wwm_fluctsplit.F90:2095:           !WRITE(STAT%FHNDL,*) IE, TRIA(IE)
wwm_fluctsplit.F90:2098:           !WRITE(STAT%FHNDL,*) IP, SI(IP)
wwm_fluctsplit.F90:2116:    write(DBG%FHNDL,*) "WARNING", __FILE__ , "Line", __LINE__
wwm_fluctsplit.F90:2117:    write(DBG%FHNDL,*) "MAXMNECON from selfe does not match self calc value. This could be problems", MAXMNECON, MNEI
wwm_fluctsplit.F90:2124:         WRITE(STAT%FHNDL,'("+TRACE......",A)') 'CALCULATE FLUCTUATION POINTER'
wwm_fluctsplit.F90:2152:           WRITE(DBG%FHNDL,*) 'COUNT_MAX=', COUNT_MAX
wwm_fluctsplit.F90:2153:           WRITE(DBG%FHNDL,*) 'MNE=', MNE
wwm_fluctsplit.F90:2218:           WRITE(STAT%FHNDL,'("+TRACE......",A)') 'SET UP SPARSE MATRIX POINTER ... COUNT NONZERO ENTRY'
wwm_fluctsplit.F90:2239:           WRITE(STAT%FHNDL,'("+TRACE......",A)') 'SET UP SPARSE MATRIX POINTER ... SETUP POINTER'
wwm_fluctsplit.F90:2420:!            WRITE(DBG%FHNDL,'(3I10,3F15.4)') IS, ID, IE, KELEM(:,IE)
wwm_fluctsplit.F90:2575:           WRITE(STAT%FHNDL,*) 'MAX. ITERATIONS USED IN ADV. SCHEME', ITER_MAX, MAXVAL(ITER_EXP)
wwm_fluctsplit.F90:2823:           WRITE(STAT%FHNDL,*) 'MAX. ITERATIONS USED IN ADV. SCHEME', ITER_MAX, MAXVAL(ITER_EXP)
wwm_fluctsplit.F90:3007:           WRITE(STAT%FHNDL,*) 'MAX. ITERATIONS USED IN ADV. SCHEME', ITER_MAX
wwm_fluctsplit.F90:3048:!                 WRITE(DBG%FHNDL,'(3I10,3F15.4)') IS, ID, IE, KELEM(:,IE)
wwm_friction.F90:45:                 WRITE(DBG%FHNDL,*) ' error in iteration fw: Madsen formulation'
wwm_friction.F90:58:           !WRITE(*,'(2I10,10F20.10)') IP, IS, KDEP, DEP(IP), CFBOT, SUM(SSBF(IS,:)), (SPSIG(IS) / SINH(MIN(20.0_rkind,KDEP)))**2 
wwm_friction.F90:97:         IF (SUM(ACLOC) .ne. SUM(ACLOC)) WRITE(*,*) SUM(ACLOC), IP
wwm_friction.F90:115:                 WRITE(DBG%FHNDL,*) ' error in iteration fw: Madsen formulation'
wwm_friction.F90:128:           !WRITE(*,'(2I10,10F20.10)') IP, IS, KDEP, DEP(IP), CFBOT, SUM(SSBF(IS,:)), (SPSIG(IS) / SINH(MIN(20.0_rkind,KDEP)))**2 
wwm_gridcf.F90:113:                      WRITE(DBG%FHNDL,*) 'WRONG ELEMENT', IE, 'WRONG NODENUMBERS', INE(:,IE)
wwm_gridcf.F90:114:                      WRITE(DBG%FHNDL,'(A40,6F15.8)') 'EDGELENGTHS OF THE WRONG ELEMENT', IEN(:,IE)
wwm_gridcf.F90:116:                     write(DBG%FHNDL,*) 'IE=', IE, ' TRIA=', TRIA(IE)
wwm_gridcf.F90:117:                     write(DBG%FHNDL,*) 'DXP1=', DXP1, ' DXP3=', DXP3
wwm_gridcf.F90:118:                     write(DBG%FHNDL,*) 'DYP1=', DYP1, ' DYP3=', DYP3
wwm_gridcf.F90:119:                     write(DBG%FHNDL,*) 'I123=', I1, I2, I3
wwm_gridcf.F90:120:                     write(DBG%FHNDL,*) 'XP,YP(I1)=', XP(I1), YP(I1)
wwm_gridcf.F90:121:                     write(DBG%FHNDL,*) 'XP,YP(I2)=', XP(I2), YP(I2)
wwm_gridcf.F90:122:                     write(DBG%FHNDL,*) 'XP,YP(I3)=', XP(I3), YP(I3)
wwm_gridcf.F90:160:                  WRITE(GRDCOR%FHNDL,'(I10)') 0
wwm_gridcf.F90:161:                  WRITE(GRDCOR%FHNDL,'(I10)') MNP
wwm_gridcf.F90:164:                      WRITE(GRDCOR%FHNDL,'(I10,3F15.8)') IP-1, XP(IP), YP(IP), DEP(IP)
wwm_gridcf.F90:168:                      WRITE(GRDCOR%FHNDL,'(I10,3F15.6)') IP-1, XP(IP), YP(IP), DEP(IP)
wwm_gridcf.F90:171:                  WRITE(GRDCOR%FHNDL,'(I10)') MNE
wwm_gridcf.F90:173:                    WRITE(GRDCOR%FHNDL,'(5I10)') INE(1,IE)-1, INE(2,IE)-1, INE(3,IE)-1, 0, IE-1
wwm_gridcf.F90:183:                  WRITE(DBG%FHNDL,*) 'The Elements in your mesh are not correctly numbered!'
wwm_gridcf.F90:184:                  WRITE(DBG%FHNDL,*) 'New mesh is written to', TRIM(GRDCOR%FNAME)
wwm_gridcf.F90:185:                  WRITE(DBG%FHNDL,*) 'There are totally', IEWRONGSUM, 'wrong Elements' 
wwm_gridcf.F90:186:                  WRITE(DBG%FHNDL,*) 'The last wrong element has the number', IEWRONG
wwm_gridcf.F90:209:                  WRITE(STAT%FHNDL,101) AVETL, TLMIN, TLMAX
wwm_gridcf.F90:211:                  WRITE(STAT%FHNDL,102) AVETA, AVETL, TLMIN, TLMAX
wwm_gridcf.F90:213:                     WRITE(STAT%FHNDL,*) ' The element area = ' 
wwm_gridcf.F90:215:                        WRITE(STAT%FHNDL,*) 'IE = ', IE, TRIA(IE)
wwm_gridcf.F90:282:         WRITE(STAT%FHNDL,*) 'RESOLUTION IN SIGMA SPACE AND FACTORS'
wwm_gridcf.F90:283:         WRITE(STAT%FHNDL,*) 'SGLOW', SGLOW
wwm_gridcf.F90:284:         WRITE(STAT%FHNDL,*) 'LOPTSIG', LOPTSIG
wwm_gridcf.F90:285:         WRITE(STAT%FHNDL,*) 'SFAC, FRINTF, FRINTH', SFAC, FRINTF, FRINTH
wwm_gridcf.F90:293:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'REL. FREQ. Distribution is =', FRINTF 
wwm_gridcf.F90:296:           WRITE(DBG%FHNDL,*) 'Freq. resolution is not optimal for Snl4'
wwm_gridcf.F90:297:           WRITE(DBG%FHNDL,'(3F15.4)') 1. + FRINTF, ABS(FRINTF - .1)/FRINTF * 100.
wwm_gridcf.F90:298:           WRITE(DBG%FHNDL,*) 'rel. freq. res. should be 1.1 is now', 1. + FRINTF, 'ERROR IS:', ABS(FRINTF - .1)/FRINTF * 100.
wwm_gridcf.F90:366:              WRITE(STAT%FHNDL,*) 'MINDIR MAXDIR', MINDIR, MAXDIR, MINDIR*RADDEG, MAXDIR*RADDEG
wwm_hotfile.F90:288:          WRITE(errmsg, *) 'Not enough data nbProc=', nbProc, ' nbMissedMode=', nbZero
wwm_hotfile.F90:442:        WRITE(HOTOUT%FHNDL) NP_TOTAL, NE_TOTAL
wwm_hotfile.F90:443:        WRITE(HOTOUT%FHNDL) MSC, MDC, FRLOW, FRHIGH
wwm_hotfile.F90:445:        WRITE(HOTOUT%FHNDL) AC2
wwm_hotfile.F90:471:            WRITE(HOTOUT%FHNDL) ACreturn
wwm_hotfile.F90:475:          WRITE(HOTOUT%FHNDL) MNP
wwm_hotfile.F90:476:          WRITE(HOTOUT%FHNDL) IPLG
wwm_hotfile.F90:477:          WRITE(HOTOUT%FHNDL) AC2
wwm_implsch2.F90:107:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) '--------- COMPUTING SOURCE TERMS ---------'
wwm_implsch2.F90:112:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'HS and TM'
wwm_implsch2.F90:113:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(10F15.7)') 4*SQRT(EMEAN(IJS)), FMEAN(IJS), SUM(FL3)
wwm_implsch2.F90:114:      !WRITE(55555,*) 4*SQRT(EMEAN(IJS)), FMEAN(IJS)
wwm_implsch2.F90:116:!      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'DIRECTIONAL PROPERTIES'
wwm_implsch2.F90:120:!          WRITE(111113,'(I10,10F15.7)') K, SPRD(IJ,K), TH(K), THWNEW(IJ) 
wwm_implsch2.F90:144:        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED BEFORE DO LOOP'
wwm_implsch2.F90:148:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER AIRSEA 1'
wwm_implsch2.F90:149:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, U10NEW(IJS), TAUW(IJS), &
wwm_implsch2.F90:162:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER SINPUT 1'
wwm_implsch2.F90:163:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL)
wwm_implsch2.F90:165:        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED'
wwm_implsch2.F90:180:        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED'
wwm_implsch2.F90:184:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER STRESSO 1'
wwm_implsch2.F90:185:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(2I10,7F15.7,I10)') IJS, IJL, SUM(FL3), &
wwm_implsch2.F90:193:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER AIRSEA 2'
wwm_implsch2.F90:194:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, U10NEW(IJS), TAUW(IJS), &
wwm_implsch2.F90:197:        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED'
wwm_implsch2.F90:204:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER SNON'
wwm_implsch2.F90:205:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL)
wwm_implsch2.F90:207:        WRITE(IU06,*) '   SUB. IMPLSCH: SNONLIN CALLED'
wwm_implsch2.F90:217:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER DISSIP' 
wwm_implsch2.F90:218:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL) 
wwm_implsch2.F90:220:        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED'
wwm_implsch2.F90:226:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER SBOTTOM' 
wwm_implsch2.F90:227:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL) 
wwm_implsch2.F90:236:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) '-------- FINISHED SOURCE TERM COMPUTATION ----------'
wwm_implsch2.F90:311:        !IF (LOUTWAM) WRITE(111113,'(4F20.10)') USNEW(IJ), FMEANWS(IJ), FMEAN(IJ)
wwm_implsch2.F90:317:          !WRITE(111113,'(4F20.10)') DELFL(M), COFRM4(M), DELT
wwm_implsch2.F90:333:      !IF (LOUTWAM) WRITE(111113,*) 'AFTER INTEGRATION' 
wwm_implsch2.F90:334:      !IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL), SUM(FL3)
wwm_implsch2.F90:396:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) '------------ INIT OF POST SOURCE TERMS --------------'
wwm_implsch2.F90:397:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(F20.10,3I10)') SUM(FL3), SIZE(FL3), IJS, IJL
wwm_implsch2.F90:414:!            WRITE(111113,'(4F20.10)') AKM1, AK2VGM1, TEMP2(IJ,M) 
wwm_implsch2.F90:448:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER MEAN PARAMETER'
wwm_implsch2.F90:449:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,5F20.10)') MIJ(IJ), AKMEAN, FMEANWS, TEMP2(IJ,MIJ(IJ)), GADIAG(IJ)
wwm_implsch2.F90:460:      !    WRITE(111113,'(I10,2F15.10)') M, FCONST(IJ,M), TEMP(IJ,M)
wwm_implsch2.F90:465:      !WRITE(111113,'(I10,3F15.10)')M,FCONST(IJ,M),TEMP(IJ,M),GADIAG(IJ)
wwm_implsch2.F90:473:      !    WRITE(111113,'(I10,10F20.10)') K, FL3(IJ,K,MIJ(IJ))
wwm_implsch2.F90:482:!            WRITE(111113,'(2I10,10F20.10)')K,M,FL3(IJ,K,M),&
wwm_implsch2.F90:489:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'BEFORE WIND INPUT 2'
wwm_implsch2.F90:490:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(8F20.10)') SUM(FL3) , SUM(FL), THWNEW(IJ)
wwm_implsch2.F90:491:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(8F20.10)') USNEW(IJ), Z0NEW(IJ), ZIDLNEW(IJ)
wwm_implsch2.F90:492:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(8F20.10)') SUM(SL), SUM(XLLWS(IJ,:,:))
wwm_implsch2.F90:506:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER SINPUT 2'
wwm_implsch2.F90:507:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL)
wwm_implsch2.F90:510:        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED AT THE END'
wwm_implsch2.F90:518:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER STRESSO 2'
wwm_implsch2.F90:519:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(2I10,7F15.8,I10)') IJS, IJL,SUM(FL3),THWNEW, USNEW, Z0NEW,ROAIRN, TAUW,SUM(SL),MIJ(IJS:IJL)
wwm_implsch2.F90:522:        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED AT THE END'
wwm_implsch2.F90:530:        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED AT THE END'
wwm_implsch2.F90:534:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER AIRSEA3'
wwm_implsch2.F90:535:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,4F15.7,I10)') &
wwm_implsch2.F90:546:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER DISSIP' 
wwm_implsch2.F90:547:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(FL3), SUM(SL) 
wwm_implsch2.F90:549:        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED AT THE END'
wwm_implsch2.F90:592:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(5F20.10)') SUM(FL3), UFRIC(IJ), Z0(IJ), CD(IJ)
wwm_implsch2.F90:593:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) '------------ END OF POST SOURCE TERMS --------------'
wwm_implsch.F90:202:      IF (LOUTWAM) WRITE(111113,*) 'HS and TM'
wwm_implsch.F90:203:      IF (LOUTWAM) WRITE(111113,'(10F15.7)') 4*SQRT(EMEAN(IJS)), FMEAN(IJS)
wwm_implsch.F90:204:      !WRITE(55555,*) 4*SQRT(EMEAN(IJS)), FMEAN(IJS)
wwm_implsch.F90:206:      IF (LOUTWAM) WRITE(111113,*) 'DIRECTIONAL PROPERTIES'
wwm_implsch.F90:210:!          WRITE(111113,'(I10,10F15.7)') K, SPRD(IJ,K), TH(K), THWNEW(IJ) 
wwm_implsch.F90:219:      IF (LOUTWAM) WRITE(111113,*) 'SOME THINKS THAT DO NOT NEED TO BE ALWAYS RECOMPUTED'
wwm_implsch.F90:220:      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJ,JU(IJS:IJL),JUMAX,MAX(NINT(XJ),1)
wwm_implsch.F90:236:        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED BEFORE DO LOOP'
wwm_implsch.F90:240:      IF (LOUTWAM) WRITE(111113,*) 'AFTER AIRSEA 1'
wwm_implsch.F90:241:      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, U10NEW(IJS), TAUW(IJS), &
wwm_implsch.F90:254:      IF (LOUTWAM) WRITE(111113,*) 'AFTER SINPUT 1'
wwm_implsch.F90:255:      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL)
wwm_implsch.F90:257:        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED'
wwm_implsch.F90:272:        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED'
wwm_implsch.F90:276:      IF (LOUTWAM) WRITE(111113,*) 'AFTER STRESSO 1'
wwm_implsch.F90:277:      IF (LOUTWAM) WRITE(111113,'(2I10,15F15.7)') IJS, IJL, SUM(FL3), &
wwm_implsch.F90:285:      IF (LOUTWAM) WRITE(111113,*) 'AFTER AIRSEA 2'
wwm_implsch.F90:286:      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, U10NEW(IJS), TAUW(IJS), &
wwm_implsch.F90:289:        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED'
wwm_implsch.F90:296:      IF (LOUTWAM) WRITE(111113,*) 'AFTER SNON'
wwm_implsch.F90:297:      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL)
wwm_implsch.F90:299:        WRITE(IU06,*) '   SUB. IMPLSCH: SNONLIN CALLED'
wwm_implsch.F90:309:      IF (LOUTWAM) WRITE(111113,*) 'AFTER DISSIP' 
wwm_implsch.F90:310:      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL) 
wwm_implsch.F90:312:        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED'
wwm_implsch.F90:318:      IF (LOUTWAM) WRITE(111113,*) 'AFTER SBOTTOM' 
wwm_implsch.F90:319:      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL) 
wwm_implsch.F90:338:        IF (LOUTWAM) WRITE(111113,'(4F20.10)') USNEW(IJ), FMEANWS(IJ), FMEAN(IJ)
wwm_implsch.F90:344:!          WRITE(111113,'(4F20.10)') DELFL(M), COFRM4(M), DELT
wwm_implsch.F90:362:!            WRITE(111113,'(4F20.10)') AKM1, AK2VGM1, TEMP2(IJ,M) 
wwm_implsch.F90:366:!      WRITE(111113,*) 'MORE TEST'
wwm_implsch.F90:377:!      WRITE(111113,'(4F20.10)')GTEMP2,FLHAB,TEMP(IJ,M)
wwm_implsch.F90:382:      IF (LOUTWAM) WRITE(111113,*) 'AFTER INTEGRATION' 
wwm_implsch.F90:383:      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL), SUM(FL3)
wwm_implsch.F90:409:        IF (LOUTWAM) WRITE(111113,*) 'AFTER MEAN PARAMETER'
wwm_implsch.F90:410:        IF (LOUTWAM) WRITE(111113,'(I10,5F20.10)') MIJ(IJ), AKMEAN, FMEANWS, TEMP2(IJ,MIJ(IJ)), GADIAG(IJ)
wwm_implsch.F90:421:      !    WRITE(111113,'(I10,2F15.10)') M, FCONST(IJ,M), TEMP(IJ,M)
wwm_implsch.F90:426:      !WRITE(111113,'(I10,3F15.10)')M,FCONST(IJ,M),TEMP(IJ,M),GADIAG(IJ)
wwm_implsch.F90:434:      !    WRITE(111113,'(I10,10F20.10)') K, FL3(IJ,K,MIJ(IJ))
wwm_implsch.F90:443:!            WRITE(111113,'(2I10,10F20.10)')K,M,FL3(IJ,K,M),&
wwm_implsch.F90:450:        IF (LOUTWAM) WRITE(111113,*) 'BEFORE WIND INPUT 2'
wwm_implsch.F90:451:        IF (LOUTWAM) WRITE(111113,'(8F20.10)') SUM(FL3) , SUM(FL), THWNEW(IJ)
wwm_implsch.F90:452:        IF (LOUTWAM) WRITE(111113,'(8F20.10)') USNEW(IJ), Z0NEW(IJ), ZIDLNEW(IJ)
wwm_implsch.F90:453:        IF (LOUTWAM) WRITE(111113,'(8F20.10)') SUM(SL), SUM(XLLWS(IJ,:,:))
wwm_implsch.F90:467:      IF (LOUTWAM) WRITE(111113,*) 'AFTER SINPUT 2'
wwm_implsch.F90:468:      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL)
wwm_implsch.F90:471:        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED AT THE END'
wwm_implsch.F90:479:      IF (LOUTWAM) WRITE(111113,*) 'AFTER STRESSO 2'
wwm_implsch.F90:480:      IF (LOUTWAM) WRITE(111113,'(I10,15F15.7)') IJS, IJL, SUM(FL3), &
wwm_implsch.F90:488:        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED AT THE END'
wwm_implsch.F90:496:        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED AT THE END'
wwm_implsch.F90:500:      IF (LOUTWAM) WRITE(111113,*) 'AFTER AIRSEA 3'
wwm_implsch.F90:501:      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, U10NEW(IJS), TAUW(IJS), &
wwm_implsch.F90:511:      IF (LOUTWAM) WRITE(111113,*) 'AFTER DISSIP' 
wwm_implsch.F90:512:      IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(FL3), SUM(SL) 
wwm_implsch.F90:514:        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED AT THE END'
wwm_inisnonlin.F90:251:        WRITE(IU06,*) '*************************************'
wwm_inisnonlin.F90:252:        WRITE(IU06,*) 'ERROR IN INISNONLIN : ICOUNT NE NINL'
wwm_inisnonlin.F90:253:        WRITE(IU06,*) 'ICOUNT= ',ICOUNT
wwm_inisnonlin.F90:254:        WRITE(IU06,*) 'NINL= ',NINL
wwm_inisnonlin.F90:255:        WRITE(IU06,*) '*************************************'
wwm_inisnonlin.F90:259:        WRITE(IU06,*) '*************************************'
wwm_inisnonlin.F90:260:        WRITE(IU06,*) 'ERROR IN INISNONLIN : IRCOUNT NE NRNL'
wwm_inisnonlin.F90:261:        WRITE(IU06,*) 'IRCOUNT= ',IRCOUNT
wwm_inisnonlin.F90:262:        WRITE(IU06,*) 'NRNL= ',NRNL
wwm_inisnonlin.F90:263:        WRITE(IU06,*) '*************************************'
wwm_initio.F90:313:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'READ SPATIAL GRID'
wwm_initio.F90:446:            write(wwmerr,*)'MAIN: nx1 wrong',i,j,nx1(i,j)
wwm_initio.F90:453:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE SETTING FHNDL'
wwm_initio.F90:457:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE READING NAMELIST'
wwm_initio.F90:475:      write(DBG%FHNDL,*) 'sum(XPtotal)=', sum(XPtotal)
wwm_initio.F90:476:      write(DBG%FHNDL,*) 'sum(YPtotal)=', sum(YPtotal)
wwm_initio.F90:514:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT SPATIAL GRID'
wwm_initio.F90:521:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT DISLIN                '
wwm_initio.F90:527:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'CHECK LOGICS                '
wwm_initio.F90:544:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET DEPTH POINTER'
wwm_initio.F90:547:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE SPECTRAL GRID'
wwm_initio.F90:551:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE BOUNDARY POINTER 1/2'
wwm_initio.F90:559:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE BOUNDARY POINTER 2/2'
wwm_initio.F90:564:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'THE FLUCTUATION SPLITTING PREPROCESSOR HAS STARTED'
wwm_initio.F90:579:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'THE FLUCTUATION SPLITTING PREPROCESSOR HAS ENDED'
wwm_initio.F90:593:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE WIND CURRENT WATERLEVEL'
wwm_initio.F90:604:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'COMPUTE THE WAVE PARAMETER'
wwm_initio.F90:611:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT ARDHUIN et al.'
wwm_initio.F90:620:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET THE INITIAL WAVE BOUNDARY CONDITION'
wwm_initio.F90:623:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET THE INITIAL CONDITION'
wwm_initio.F90:626:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET BOUNDARY CONDITIONS'
wwm_initio.F90:629:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT STATION OUTPUT'
wwm_initio.F90:632:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'WRITING INITIAL TIME STEP'
wwm_initio.F90:647:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'OPEN PIPES FOR COUPLING'
wwm_initio.F90:670:        WRITE(DBG%FHNDL,*) 'MSC_SELFE', MSC_SELFE
wwm_initio.F90:671:        WRITE(DBG%FHNDL,*) 'MSC', MSC
wwm_initio.F90:672:        WRITE(DBG%FHNDL,*) 'MDC_SELFE', MDC_SELFE
wwm_initio.F90:673:        WRITE(DBG%FHNDL,*) 'MDC', MDC
wwm_initio.F90:679:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE_WWM'
wwm_initio.F90:681:      WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'CPU Time for the preprocessing', TIME2-TIME1
wwm_initio.F90:923:         WRITE(STAT%FHNDL,*) 'START WAVE PARAMETER'
wwm_initio.F90:926:         WRITE(STAT%FHNDL,*) 'GRADDEP'
wwm_initio.F90:929:         WRITE(STAT%FHNDL,*) 'GRADCURT'
wwm_initio.F90:932:         WRITE(STAT%FHNDL,*) 'BASIC'
wwm_initio.F90:936:         WRITE(STAT%FHNDL,*) 'WAVEKCG'
wwm_initio.F90:952:           WRITE(STAT%FHNDL,'("+TRACE...",A)')'COMPUTING NONLINEAR COEFFICIENTS' 
wwm_initio.F90:954:           WRITE(STAT%FHNDL,'("+TRACE...",A)')'COMPUTING NONLINEAR COEFFICIENTS'
wwm_initio.F90:959:             WRITE(STAT%FHNDL,'("+TRACE...",A)')'READING STRESS TABLES'
wwm_initio.F90:973:             WRITE(STAT%FHNDL,'("+TRACE...",A)')'COMPUTING STRESS TABLES'
wwm_initio.F90:975:             WRITE(STAT%FHNDL,'("+TRACE...",A)')'COMPUTING HF TABLES'
wwm_initio.F90:979:           WRITE(STAT%FHNDL,'("+TRACE...",A)')'INITIALIZING STRESS ARRAYS'
wwm_initio.F90:1023:              WRITE(2001) 1.
wwm_initio.F90:1024:              WRITE(2001) (TMPPAR(3,IP), TMPPAR(2,IP), TMPPAR(1,IP), IP = 1, MNP)
wwm_initio.F90:1139:         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
wwm_initio.F90:1143:         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
wwm_initio.F90:1147:         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
wwm_initio.F90:1152:         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
wwm_initio.F90:1156:         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
wwm_initio.F90:1160:         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
wwm_initio.F90:1166:         WRITE(DBG%FHNDL, *) 'THR=', THR
wwm_initio.F90:1167:         WRITE(DBG%FHNDL, *) 'THR8=', THR8
wwm_initio.F90:1172:         WRITE(STAT%FHNDL,*) 'Input Filename   =', TRIM(INP%FNAME)
wwm_initio.F90:1173:         WRITE(STAT%FHNDL,*) 'Check Filename   =', TRIM(CHK%FNAME)
wwm_initio.F90:1174:         WRITE(STAT%FHNDL,*) 'Qstea Filename   =', TRIM(QSTEA%FNAME)
wwm_initio.F90:1175:         WRITE(STAT%FHNDL,*) 'Iobp Filename    =', TRIM(IOBPOUT%FNAME)
wwm_initio.F90:1176:         WRITE(STAT%FHNDL,*) 'Iobpd Filename   =', TRIM(IOBPDOUT%FNAME)
wwm_initio.F90:1177:         WRITE(STAT%FHNDL,*) 'WindDbg Filename =', TRIM(WINDBG%FNAME)
wwm_initio.F90:1249:        WRITE(DBG%FHNDL,*) 'SEARCHING FOR STATION ACROSS RANKS', myrank
wwm_initio.F90:1258:            WRITE(DBG%FHNDL,'(A10,I10,A20,I10,A15,2I10)') 'MYRANK', MYRANK, 'STATION =',I, 'IN ELEMENT =', IELG(STATION(I)%ELEMENT), STATION(I)%IFOUND
wwm_initio.F90:1265:            WRITE(DBG%FHNDL,'(A10,I10,A20,I10,A15,2I10)') 'MYRANK', MYRANK, 'STATION =',I, 'IN ELEMENT =', STATION(I)%ELEMENT, STATION(I)%IFOUND
wwm_initio.F90:1274:            WRITE(DBG%FHNDL,'(A30,3I10)') 'SUM OF THE FOUND STATIONS MYRANK', MYRANK, I, STATION(I)%ISUM
wwm_initio.F90:1289:              WRITE(DBG%FHNDL,'(A20,I10,A10,2F15.8)') 'STATION NOT FOUND', I, STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD
wwm_initio.F90:1291:              WRITE(DBG%FHNDL,'(A25,I10,A10,2F15.8)') 'STATION FOUND    ', I, STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD
wwm_initio.F90:1312:          WRITE(STAT%FHNDL,*) 'FINDING ELEMENT CONNECTED TO STATION'
wwm_initio.F90:1317:              WRITE(STAT%FHNDL,*) STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD, ' is out of mesh !'
wwm_initio.F90:1324:              WRITE(STAT%FHNDL,*)'Site    ',STATION(I)%NAME, STATION(I)%XCOORD,STATION(I)%YCOORD,STATION(I)%IFOUND
wwm_initio.F90:1338:            WRITE(DBG%FHNDL,*) 'CUT-OFF FREQ. OF STATION =', IP, STATION(IP)%CUTOFF, 'RAD - IS =', STATION(IP)%ISMAX
wwm_initio.F90:1342:      WRITE(STAT%FHNDL,'("+TRACE...",A)')'FINISHED WITH INIT_STATION_OUTPUT'
wwm_initio.F90:1519:        !WRITE(*,*) WBMSC, WBMDC
wwm_initio.F90:1535:          !WRITE(*,*) SPEG(IS,1,1), SDIR(IS,1), SPRD(IS,1)
wwm_initio.F90:1567:        WRITE(STAT%FHNDL,*) 'START READSPEC2D_WW3_INIT_SPEC'
wwm_initio.F90:1568:        WRITE(STAT%FHNDL,*) 'LABEL, MSC_WW3,MDC_WW3, NP_WW3, GNAME'
wwm_initio.F90:1569:        WRITE(STAT%FHNDL,*) LABEL, MSC_WW3,MDC_WW3, NP_WW3, GNAME
wwm_initio.F90:1571:        WRITE(STAT%FHNDL,*)'DIRECTION NUMBER IN WW3 SPECTRUM:',MDC_WW3
wwm_initio.F90:1572:        WRITE(STAT%FHNDL,*)'FREQUENCY NUMBER IN WW3 SPECTRUM:',MSC_WW3
wwm_initio.F90:1573:        WRITE(STAT%FHNDL,'("+TRACE...",A)')'DONE READSPEC2D_WW3_INIT_SPEC'
wwm_initio.F90:1601:        WRITE(STAT%FHNDL,*)'START READSPEC2D_WW3_INIT_TIME'
wwm_initio.F90:1623:            WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 1, WHICH IS NICE'
wwm_initio.F90:1628:!            WRITE(STAT%FHNDL,'(A10,7F15.4)') PID,TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7
wwm_initio.F90:1632:              WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 2, WHICH IS NOT GOOD'
wwm_initio.F90:1639:              WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 3, WHICH IS NOT GOOD'  
wwm_initio.F90:1647:        WRITE(STAT%FHNDL,*) 'NUMBER OF BUOYS', NP_WW3
wwm_initio.F90:1648:        WRITE(STAT%FHNDL,*) 'NUMBER OF TIME STEPS IN FILE', MAXSTEP_WW3 
wwm_initio.F90:1672:          WRITE(CTIME1,*) ITIME(IT,1) 
wwm_initio.F90:1676:          WRITE(STAT%FHNDL,*) IT, TIMESTRING, BND_TIME_ALL_FILES(1,IT)
wwm_initio.F90:1685:        WRITE(STAT%FHNDL,*)'MIN. FREQ. IN WW3 SPECTRUM:',FQ_WW3(1)
wwm_initio.F90:1686:        WRITE(STAT%FHNDL,*)'MAX. FREQ. IN WW3 SPECTRUM:',FQ_WW3(MSC_WW3)
wwm_initio.F90:1687:        WRITE(STAT%FHNDL,*)'NUMBER OF TIME STEPS',MAXSTEP_WW3
wwm_initio.F90:1688:        WRITE(STAT%FHNDL,*)'TIME INCREMENT IN SPECTRAL FILE', DTBOUND_WW3
wwm_initio.F90:1689:        WRITE(STAT%FHNDL,*)'FIRST TIME STEP IN WW3 SPECTRUM FILE:',BND_TIME_ALL_FILES(1,1)
wwm_initio.F90:1690:        WRITE(STAT%FHNDL,*)'BEGING TIME, END TIME and DELT of wave boundary', SEBO%BMJD, SEBO%EMJD, SEBO%DELT
wwm_initio.F90:1691:        WRITE(STAT%FHNDL,*)'BEGING TIME, END TIME and DELT of simulation', MAIN%BMJD, MAIN%EMJD, MAIN%DELT
wwm_initio.F90:1692:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE READSPEC2D_WW3INIT2'
wwm_initio.F90:1702:        write( string, '(I6)' ) INTIN 
wwm_initio.F90:1757:!            write(*,*) sum(SPECOUT_SGLE)
wwm_initio.F90:1766:!                write(*,*) 4*sqrt(m0), m1, m2, df
wwm_initio.F90:1783:      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE READSPEC2D_WW3'
wwm_initio.F90:1821:        !  WRITE(STAT%FHNDL,*) 'ORIG WW3 SUM SPEC', IBWW3, SUM(SPEC_WW3_UNSORT(:,:,IBWW3))
wwm_initio.F90:1837:!          WRITE(STAT%FHNDL,*) 'AFTER SORTING', IBWW3, SUM(SPEC_WW3(:,:,IBWW3))
wwm_initio.F90:1844:!          WRITE(STAT%FHNDL,*)'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
wwm_initio.F90:1845:!          WRITE(STAT%FHNDL,*)'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
wwm_initio.F90:1846:!          WRITE(STAT%FHNDL,*)'WW3 spectra does not encompass the whole WWM spectra, please carefully check if this makes sense for your simulations'
wwm_initio.F90:1849:!          WRITE(STAT%FHNDL,*)'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
wwm_initio.F90:1850:!          WRITE(STAT%FHNDL,*)'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
wwm_initio.F90:1855:          WRITE(DBG%FHNDL,*) 'SUMS AFTER INTERPOLATION', SUM(SPEC_WW3),SUM(SPEC_WWM)
wwm_initio.F90:1874:                !WRITE(STAT%FHNDL,*)'XP_WWM =',XP_WWM,'XP_WW3 =',XP_WW3(IBWW3)
wwm_initio.F90:1875:                !WRITE(STAT%FHNDL,*)'YP_WWM =',YP_WWM,'YP_WW3 =',YP_WW3(IBWW3)
wwm_initio.F90:1878:                !WRITE(STAT%FHNDL,*) 'orig', IBWW3, INDBWW3(IBWW3), DIST(IBWW3)
wwm_initio.F90:1882:              !  WRITE(STAT%FHNDL,*) 'sorted', IBWW3, INDBWW3(IBWW3), DIST(IBWW3)
wwm_initio.F90:1885:              !WRITE(STAT%FHNDL,'(A20, 2F20.5,3F30.10)') ' AFTER INTERPOLATION ', INDBWW3(1), INDBWW3(2), sum(SPEC_WWM(:,:,INT(INDBWW3(1)))), sum(SPEC_WWM(:,:,INT(INDBWW3(2)))), SUM(WBACOUT(:,:,IB))
wwm_initio.F90:1897:        !  WRITE(STAT%FHNDL,*) 'SUM OF WBAC', IB, SUM(WBACOUT(:,:,IB)) 
wwm_initio.F90:1900:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE GETWW3SPECTRA'
wwm_initio.F90:1912:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START WITH INIT_BINARY_WW3_SPECTRA'
wwm_initio.F90:1917:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE WITH INIT_BINARY_WW3_SPECTRA'
wwm_initio.F90:1970:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING SPECTRALINT'
wwm_initio.F90:1978:          WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'BEFORE INTERPOLATION', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM) 
wwm_initio.F90:2010:          WRITE(STAT%FHNDL,*) 'POINT NUMBER', IP
wwm_initio.F90:2011:          WRITE(STAT%FHNDL,'(A10,2F20.10,A10,2F20.10)') 'M1 = ',M1_WW3, M1_WWM, 'M2 = ',M2_WW3, M2_WWM
wwm_initio.F90:2015:          WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'AFTER INTERPOLATION', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM)
wwm_initio.F90:2026:        WRITE(STAT%FHNDL,*)'CHECKING INTEGRATED PARAMETERS AFTER JACOBIAN'
wwm_initio.F90:2050:          WRITE(STAT%FHNDL,*) 'POINT NUMBER', IP
wwm_initio.F90:2051:          WRITE(STAT%FHNDL,'(A10,2F20.10,A10,2F20.10)') 'M1 = ',M1_WW3, M1_WWM, 'M2 = ',M2_WW3, M2_WWM
wwm_initio.F90:2055:          WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'AFTER JACOBIAN', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM)
wwm_initio.F90:2221:              WRITE(STAT%FHNDL,*) 'INITIALIZING WAVE BOUNDARY IT =', IT
wwm_initio.F90:2222:              WRITE(STAT%FHNDL,*) 'SUM OF WAVE ACTION', SUM(WBACOLD), SUM(WBAC) 
wwm_initio.F90:2241:              WRITE(STAT%FHNDL,*) IFILE, IT, SUM(NDT_BND_FILE(1:IFILE-1)), NINT(DTMP/SEBO%DELT), SEBO%DELT
wwm_input.F90:7:# define wwm_print_namelist(xinp) IF (myrank.eq.0) WRITE(CHK%FHNDL, NML=xinp)
wwm_input.F90:9:# define wwm_print_namelist(xinp) WRITE(CHK%FHNDL, NML=xinp)
wwm_input.F90:486:           WRITE(DBG%FHNDL,*) 'OUT_STATION%FNAME=', TRIM(OUT_STATION%FNAME)
wwm_input.F90:487:           WRITE(DBG%FHNDL,*) 'OUT_HISTORY%FNAME=', TRIM(OUT_HISTORY%FNAME)
wwm_input.F90:578:                WRITE(wwmerr,*) 'Station ', I, ' has incorrect name'
wwm_input.F90:589:           WRITE(DBG%FHNDL,*) 'STATION X and Y Coordinates'
wwm_input.F90:590:           WRITE(DBG%FHNDL,*) STATION_P%XCOORD
wwm_input.F90:591:           WRITE(DBG%FHNDL,*) STATION_P%YCOORD
wwm_input.F90:592:           WRITE(DBG%FHNDL,*) 'STATION Names'
wwm_input.F90:593:           WRITE(DBG%FHNDL,*) STATION_P%NAME
wwm_input.F90:623:           WRITE(DBG%FHNDL,*) 'STATION X and Y Coordinates'
wwm_input.F90:624:           WRITE(DBG%FHNDL,*) STATION_P%XCOORD
wwm_input.F90:625:           WRITE(DBG%FHNDL,*) STATION_P%YCOORD
wwm_input.F90:626:           WRITE(DBG%FHNDL,*) 'STATION Names'
wwm_input.F90:627:           WRITE(DBG%FHNDL,*) STATION_P%NAME
wwm_input.F90:712:             WRITE(DBG%FHNDL) LSPHE, ICS
wwm_input.F90:718:             WRITE(DBG%FHNDL) LSPHE, ICS
wwm_input.F90:746:           WRITE(wwmerr,*)'You are running in less than 1d or LCPL = .F.',&
wwm_input.F90:801:           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'HOTFILE is used as Initital Condition'
wwm_input.F90:823:         WRITE(STAT%FHNDL,'("+TRACE...",A10,I5)') 'BOUNDARY FILE FORMAT IS', IBOUNDFORMAT
wwm_input.F90:1425:             WRITE(DBG%FHNDL,*) 'MAIN%DTCOUP=', MAIN%DTCOUP
wwm_input.F90:1426:             WRITE(DBG%FHNDL,*) 'MAIN%DELT=', MAIN%DELT
wwm_input.F90:1427:             WRITE(DBG%FHNDL,*) 'TEST=', TEST
wwm_input.F90:1428:             WRITE(DBG%FHNDL,*) 'TIME STEP OF THE WAVEMODELL CANNOT BE DiVIDIED WITHOUT A REST'
wwm_input.F90:1437:         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SWTICHES FOR THE LIMTER'
wwm_input.F90:1438:         WRITE(STAT%FHNDL,*) 'LLIMT', LLIMT
wwm_input.F90:1439:         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ACTIVATED SOURCE TERMS'
wwm_input.F90:1440:         WRITE(STAT%FHNDL,*) 'MESIN', MESIN
wwm_input.F90:1441:         WRITE(STAT%FHNDL,*) 'MESNL', MESNL
wwm_input.F90:1442:         WRITE(STAT%FHNDL,*) 'MESBR', MESBR
wwm_input.F90:1443:         WRITE(STAT%FHNDL,*) 'MESDS', MESDS
wwm_input.F90:1444:         WRITE(STAT%FHNDL,*) 'MESTR', MESTR
wwm_input.F90:1447:           WRITE(DBG%FHNDL,*) 'YOU MUST USE EITHER UNSTEADY OR STEADY WIND'
wwm_input.F90:1448:           WRITE(DBG%FHNDL,*) 'PLEASE CHECK CODE EXITS'
wwm_input.F90:1453:           WRITE(DBG%FHNDL,*) 'YOU MUST USE EITHER UNSTEADY OR STEADY CURRENTS'
wwm_input.F90:1454:           WRITE(DBG%FHNDL,*) 'PLEASE CHECK CODE EXITS'
wwm_input.F90:1459:           WRITE(DBG%FHNDL,*) 'YOU MUST USE EITHER UNSTEADY OR STEADY CURRENTS'
wwm_input.F90:1460:           WRITE(DBG%FHNDL,*) 'PLEASE CHECK CODE EXITS'
wwm_input.F90:1465:           WRITE(DBG%FHNDL,*) 'YOU MUST USE EITHER UNSTEADY OR STEADY CURRENTS'
wwm_input.F90:1466:           WRITE(DBG%FHNDL,*) 'PLEASE CHECK CODE EXITS'
wwm_input.F90:1550:        WRITE(STAT%FHNDL,*) 'Serial current Condition -----------'
wwm_input.F90:1551:        WRITE(STAT%FHNDL,*) SECU%BEGT, SECU%ENDT, SECU%ISTP, SECU%TOTL/3600.0, SECU%DELT
wwm_input.F90:1612:        WRITE(STAT%FHNDL,*) 'Serial water level Condition -----------'
wwm_input.F90:1613:        WRITE(STAT%FHNDL,*) SEWL%BEGT, SEWL%ENDT, SEWL%ISTP, SEWL%TOTL/3600.0, SEWL%DELT
wwm_input.F90:1642:        WRITE(STAT%FHNDL,*) WAV%FHNDL, WAV%FNAME, BND%FHNDL, BND%FNAME
wwm_input.F90:1651:        WRITE(STAT%FHNDL,*) NUM_NETCDF_FILES_BND
wwm_input.F90:1653:        WRITE(STAT%FHNDL,*) 'NUM_NETCDF_FILES_BND', NUM_NETCDF_FILES_BND
wwm_input.F90:1671:          WRITE(STAT%FHNDL,'(I10,10X,5A30)') IFILE, NETCDF_FILE_NAMES_BND(IFILE,:)
wwm_input.F90:1677:          write(STAT%FHNDL,*) ifile, TRIM(NETCDF_FILE_NAMES_BND(IFILE,1))
wwm_input.F90:1691:          write(STAT%FHNDL,*) IFILE, NDT_BND_FILE(IFILE)
wwm_input.F90:1714:        WRITE(STAT%FHNDL,*) 'Number of Gridpoints', NDX_BND, NDY_BND
wwm_input.F90:1745:        write(STAT%FHNDL,*) NUM_NETCDF_FILES_BND
wwm_input.F90:1748:          write(STAT%FHNDL,*) it, NDT_BND_FILE(it)
wwm_input.F90:1751:        WRITE(STAT%FHNDL,*) NDT_BND_ALL_FILES, NDT_BND_FILE
wwm_input.F90:1777:!             WRITE(*,*) '19000101.000000', DTMP1, chrdate
wwm_input.F90:1779:!             WRITE(*,*) '19900101.000000', DTMP2, chrdate
wwm_input.F90:1781:!             WRITE(*,*) '00000000.000000', 0.0_rkind, chrdate
wwm_input.F90:1782:!             WRITE(*,*) BND_TIME_ALL_FILES(1,1), DT_DIFF_19901900
wwm_input.F90:1783:!             IF (IT == 1 .AND. IFILE ==1) WRITE(*,*) DTMP1, DTMP2, DTMP1+DT_DIFF_19901900
wwm_input.F90:1784:!             IF (IT == 1 .AND. IFILE ==1) WRITE(*,*) IFILE, IT, BND_TIME(IT), chrdate
wwm_input.F90:1790:        write(STAT%FHNDL,*) SEBO%DELT, BND_TIME_ALL_FILES(1,2), BND_TIME_ALL_FILES(1,1)
wwm_input.F90:1841:         WRITE(DBG%FHNDL,*) IT, IFILE, 'READING GLOBAL DATA'
wwm_input.F90:1976:           WRITE(3012) TIME
wwm_input.F90:1977:           WRITE(3012) (U(IP), V(IP), H(IP), IP = 1, NDX_BND*NDY_BND)
wwm_input.F90:2048:            !WRITE(*,*) 'XP YP', XP(IWBNDLC(IP)), YP(IWBNDLC(IP)), IWBNDLC(IP)
wwm_input.F90:2049:            !WRITE(*,*) LEN_X, LEN_Y, OFFSET_X, OFFSET_Y, XP(IWBNDLC(IP)), YP(IWBNDLC(IP))
wwm_input.F90:2059:            !WRITE(*,*) WX1, WX2, WX3, WX4
wwm_input.F90:2246:            !WRITE(*,'(I10,8F15.4)') IP, VAL(:,IP)
wwm_input.F90:2253:           WRITE(4013) RTIME
wwm_input.F90:2254:           WRITE(4013) (VAL(3,IP), VAL(2,IP), VAL(4,IP), IP = 1, MNP)
wwm_main.F90:50:         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING WWM_II'
wwm_main.F90:58:           WRITE(DBG%FHNDL,*) ' STARTING WWM FROM SELFE ',  SUM(AC2)
wwm_main.F90:78:           WRITE(DBG%FHNDL,*) 'MAIN%DELT=', MAIN%DELT, ' in wwminput.nml'
wwm_main.F90:79:           WRITE(DBG%FHNDL,*) 'But nstep_wwm*dt=', DT_PROVIDED
wwm_main.F90:80:           WRITE(DBG%FHNDL,*) 'nstep_wwm=', NSTEPWWM
wwm_main.F90:81:           WRITE(DBG%FHNDL,*) '       dt=', DT_SELFE
wwm_main.F90:86:           WRITE(DBG%FHNDL,*) ' FIRST SUM IN MAIN ',  SUM(AC2)
wwm_main.F90:90:         WRITE(STAT%FHNDL,'("+TRACE...",A)') ' ---- ALL CHECKS DONE'
wwm_main.F90:213:           WRITE(DBG%FHNDL,*) ' AFTER SETTING BOUNDARY CONDITION IN MAIN ',  SUM(AC2)
wwm_main.F90:227:           WRITE(DBG%FHNDL,*) ' BEFORE COMPUTE ',  SUM(AC2)
wwm_main.F90:231:         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE'
wwm_main.F90:242:           WRITE(DBG%FHNDL,*) ' AFTER COMPUTE ',  SUM(AC2)
wwm_main.F90:250:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'FINISHED COMPUTE nth call to WWM', SIMUTIME
wwm_main.F90:299:           WRITE(DBG%FHNDL,*) ' END OF MAIN ',  SUM(AC2)
wwm_main.F90:304:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----TOTAL TIMINGS-----'
wwm_main.F90:305:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPARATION        ', TIME2-TIME1
wwm_main.F90:306:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'INTEGRATION        ', TIME3-TIME2
wwm_main.F90:307:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'OUTPUT TO SELFE    ', TIME4-TIME3
wwm_main.F90:308:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'RADIATION STRESSES ', TIME5-TIME4
wwm_main.F90:309:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'NAN CHECK          ', TIME6-TIME5
wwm_main.F90:310:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'TOTAL TIME         ', TIME6-TIME1
wwm_main.F90:311:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '------END-TIMINGS-  ---'
wwm_main.F90:313:         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'FINISHED WITH WWM', SIMUTIME
wwm_main.F90:327:          WRITE(DBG%FHNDL,*) 'NaN in WINDX', IP, WINDXY(IP,1) 
wwm_main.F90:331:          WRITE(DBG%FHNDL,*) 'NaN in WINDY', IP, WINDXY(IP,2) 
wwm_main.F90:335:          WRITE(DBG%FHNDL,*) 'NaN in WATLEV', IP, WATLEV(IP) 
wwm_main.F90:339:          WRITE(DBG%FHNDL,*) 'NaN in WATLEV', IP, WATLEV(IP)
wwm_main.F90:343:          WRITE(DBG%FHNDL,*) 'NaN in CURTX', IP, CURTXY(IP,1)
wwm_main.F90:347:          WRITE(DBG%FHNDL,*) 'NaN in CURTY', IP, CURTXY(IP,2)
wwm_main.F90:362:            WRITE(DBG%FHNDL,*) 'NaN in OUTT_INTPAR', IP, I, OUTT_INTPAR(IP,I)
wwm_main.F90:368:            WRITE(DBG%FHNDL,*) 'NaN in WIND_INTPAR', IP, I, WIND_INTPAR(IP,I)
wwm_main.F90:374:            WRITE(DBG%FHNDL,*) 'NaN in WWAVE_FORCE', IP, I, WWAVE_FORCE(I,IP,1), WWAVE_FORCE(I,IP,2)
wwm_main.F90:398:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING NON_STEADY' 
wwm_main.F90:436:      IF (myrank.eq.0) WRITE(*,101)  K, MAIN%ISTP, RTIME
wwm_main.F90:438:      WRITE(*,101)  K, MAIN%ISTP, RTIME
wwm_main.F90:463:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6,A20)') '-----SIMULATION TIME-----        ', MAIN%TMJD, CTIME
wwm_main.F90:464:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----TOTAL RUN TIMES-----        '
wwm_main.F90:465:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPROCESSING                    ', TIME2-TIME1
wwm_main.F90:466:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'INTEGRATION                      ', TIME3-TIME2
wwm_main.F90:467:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'WAVE SETUP                       ', TIME4-TIME3
wwm_main.F90:468:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'POSTPROCESSING                   ', TIME5-TIME4
wwm_main.F90:469:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CHECK STEADY                     ', TIME6-TIME5
wwm_main.F90:470:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'TOTAL TIME                       ', TIME6-TIME1
wwm_main.F90:471:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
wwm_main.F90:478:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'LEAVING NON_STEADY'
wwm_main.F90:533:            WRITE(QSTEA%FHNDL,'(3I10,5F15.8)') K, IT, NQSITER, CONV1, CONV2, CONV3, CONV4, CONV5
wwm_main.F90:535:            if (myrank == 0) WRITE(QSTEA%FHNDL,'(3I10,5F15.8)') K, IT, NQSITER, CONV1, CONV2, CONV3, CONV4, CONV5
wwm_main.F90:547:      IF (myrank == 0) WRITE(STAT%FHNDL,101)  K, MAIN%ISTP, RTIME*DAY2SEC
wwm_main.F90:551:      WRITE(STAT%FHNDL,101)  K, MAIN%ISTP, RTIME*DAY2SEC
wwm_main.F90:652:             WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'WRITING HOTFILE INTERNAL TIME', RTIME
wwm_main.F90:829:      WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----TOTAL TIME IN PROG-----', TIME2-TIME1
wwm_m_fileio.F90:178:  if(i_print >0) write(i_out,*) 'Z_FILEIO: Incorrect value for IUFIND:',iufind
wwm_m_fileio.F90:186:if(i_print>=1) write(i_out,*) 'Z_FILEIO/A:',trim(filename),' ',qual,iunit,iostat
wwm_m_fileio.F90:189:  if(i_print > 0) write(i_out,*) 'Incorrect file qualifier'
wwm_m_fileio.F90:206:  if(i_print >=2) write(i_out,*) 'Z_FILEIO  file exists?:',trim(filename),':',lexist
wwm_m_fileio.F90:232:    if(i_print >=2) write(i_out,*) 'Z_FILEIO: File exists:',trim(filename)
wwm_m_fileio.F90:235:      if(i_print >=2) write(i_out,*) 'Z_FILEIO: File is opened:',trim(filename)
wwm_m_fileio.F90:241:      if(i_print >=2) write(i_out,*) 'Z_FILEIO: File is connected to unit:', junit
wwm_m_fileio.F90:247:      if(i_print >=2) write(i_out,*) 'Z_FILEIO: File is not connected to a unit number'
wwm_m_fileio.F90:249:        if(i_print >=2) write(i_out,*) 'Z_FILEIO: Assign user defined unit number:',iunit
wwm_m_fileio.F90:252:        if(i_print >=2) write(i_out,*) 'Z_FILEIO: New unit number IUNIT:',iunit
wwm_m_fileio.F90:268:       write(i_out,*) 'Z_FILEIO: File does not exist !'
wwm_m_fileio.F90:269:       write(i_out,*) 'Z_FILEIO: Qual:',qual(1:1)
wwm_m_fileio.F90:275:        if(i_print >=1) write(i_out,*) 'Z_FILEIO: New unit number IUNIT:',iunit
wwm_m_fileio.F90:297:      if(i_print>=2) write(i_out,*) 'Z_FILEIO: File cannot be opened because it does not exist'
wwm_m_fileio.F90:305:if(i_print>=1) write(i_out,*) 'Z_FILEIO/Z:',trim(filename),' ',qual,iunit,iostat
wwm_m_fileio.F90:441:  write(i_out,*) 'Z_FLUNIT: forbidden     :',lu_not
wwm_m_fileio.F90:442:  write(i_out,*) 'Z_FLUNIT: lu_min lu_max :',lu_min,lu_max
wwm_m_fileio.F90:449:  write(i_out,*) 'Z_FLUNIT: Incorrect boundaries for LU_MIN & LU_MAX:',&
wwm_m_fileio.F90:469:       if(i_print >= 1) write(i_out,*) 'Z_FLUNIT: a forbidden unit number was encountered:',junit
wwm_m_fileio.F90:482:  write(i_out,*) 'ERROR in Z_FLUNIT: No free unit number could be found'
wwm_mjdv2.F90:52:         WRITE(STIME(1:4),'(I4.4)') IY
wwm_mjdv2.F90:53:         WRITE(STIME(5:6),'(I2.2)') IM
wwm_mjdv2.F90:54:         WRITE(STIME(7:8),'(I2.2)') ID
wwm_mjdv2.F90:55:         WRITE(STIME(9:9),'(A)') '.'
wwm_mjdv2.F90:56:         WRITE(STIME(10:11),'(I2.2)') IH
wwm_mjdv2.F90:57:         WRITE(STIME(12:13),'(I2.2)') IMIN
wwm_mjdv2.F90:58:         WRITE(STIME(14:15),'(I2.2)') ISEC
wwm_mjdv2.F90:76:            WRITE(DBG%FHNDL,*) 'ERROR WRONG UNIT, UNIT = ', UNITT
wwm_mod_xnl4v5.F90:639:   write(iscreen,'(a,i4,2f10.2)') 'XNL_INIT iq_grid dstep dgap:',iq_grid,dstep,dgap
wwm_mod_xnl4v5.F90:646:   write(iscreen,'(a,i4,2f10.2)') 'XNL_INIT iq_grid dird(1) dird(n):',iq_grid,dird(1),dird(ndir)
wwm_mod_xnl4v5.F90:695:write(luq_log,'(2a,i4)') 'XNL_INIT: ',trim(qbase)//'.log connected to :',luq_log
wwm_mod_xnl4v5.F90:696:write(luq_log,'(2a,i4)') 'XNL_INIT: ',trim(qbase)//'.prt connected to :',luq_prt
wwm_mod_xnl4v5.F90:697:write(luq_log,'(2a,i4)') 'XNL_INIT: ',trim(qbase)//'.tst connected to :',luq_tst
wwm_mod_xnl4v5.F90:699:write(luq_log,'(2a,i4)') 'XNL_INIT: ',trim(qbase)//'.int connected to :',luq_int
wwm_mod_xnl4v5.F90:700:write(luq_log,'(2a,i4)') 'XNL_INIT: ',trim(qbase)//'.trf connected to :',luq_trf
wwm_mod_xnl4v5.F90:701:write(luq_log,'(2a,i4)') 'XNL_INIT: ',trim(qbase)//'.t13 connected to :',luq_t13
wwm_mod_xnl4v5.F90:703:write(luq_prt,'(a)') '---------------------------------------------------------------'
wwm_mod_xnl4v5.F90:704:write(luq_prt,'(a)') trim(q_version)
wwm_mod_xnl4v5.F90:705:write(luq_prt,'(a)') 'Solution of Boltzmann integral using Webb/Resio/Tracy method'
wwm_mod_xnl4v5.F90:706:write(luq_prt,'(a)') '---------------------------------------------------------------'
wwm_mod_xnl4v5.F90:707:write(luq_prt,*)
wwm_mod_xnl4v5.F90:708:write(luq_prt,'(a)') 'Initialisation'
wwm_mod_xnl4v5.F90:709:write(luq_prt,*)
wwm_mod_xnl4v5.F90:711:if(iproc >=0) write(luq_prt,'(a,i5)') '(MPI) processor number:',iproc
wwm_mod_xnl4v5.F90:753:    write(luq_err,'(a,e12.5,f10.2)') 'Incorrect depth & minimum:',q_depth,q_mindepth
wwm_mod_xnl4v5.F90:764:    write(luq_prt,'(a)') 'XNL_INIT: For deep water only one grid suffices'
wwm_mod_xnl4v5.F90:903:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:904:  write(luq_prt,'(a,i4,f16.3,i4)') 'XNL_MAIN: Input arguments: iquad depth iproc:',&
wwm_mod_xnl4v5.F90:958:    if(iq_prt >=1) write(luq_prt,'(a,f7.4)') 'XNL_MAIN depth scale factor:',q_dfac
wwm_mod_xnl4v5.F90:967:    write(luq_prt,'(a)')        'XNL_MAIN: Conservation checks'
wwm_mod_xnl4v5.F90:968:    write(luq_prt,'(a,4e13.5)') 'XNL_MAIN: E/A/MOMX/MOMY:',sum_e,sum_a,sum_mx,sum_my
wwm_mod_xnl4v5.F90:976:  write(luq_log,*)
wwm_mod_xnl4v5.F90:977:  write(luq_log,'(a,i4)') 'XNL_MAIN: Number of warnings:',iq_warn
wwm_mod_xnl4v5.F90:978:  write(luq_log,'(a,i4)') 'XNL_MAIN: Number of errors  :',iq_err
wwm_mod_xnl4v5.F90:1073:  write(iscreen,'(a,2i4,2f8.2)') 'Q_ADDTAIL nsig,isig_tail,qf_tail:',nsig,isig_tail,pf_tail
wwm_mod_xnl4v5.F90:1074:  write(iscreen,'(a,3f8.2)') 'Q_ADDTAIL qf_tail, x_tail d_tail:',pf_tail,x_tail,d_tail
wwm_mod_xnl4v5.F90:1176:if(iq_log >=1) write(luq_log,'(a,4i4)') &
wwm_mod_xnl4v5.F90:1320:  write(luq_log,'(a)')    'Q_ALLOCATE: size of arrays'
wwm_mod_xnl4v5.F90:1321:  write(luq_log,'(a,i4)') 'Q_ALLOCATE: mkq  :',mkq
wwm_mod_xnl4v5.F90:1322:  write(luq_log,'(a,i4)') 'Q_ALLOCATE: maq  :',maq
wwm_mod_xnl4v5.F90:1323:  write(luq_log,'(a,i4)') 'Q_ALLOCATE: nkq  :',nkq
wwm_mod_xnl4v5.F90:1324:  write(luq_log,'(a,i4)') 'Q_ALLOCATE: naq  :',naq
wwm_mod_xnl4v5.F90:1325:  write(luq_log,'(a,i4)') 'Q_ALLOCATE: mlocus:',mlocus
wwm_mod_xnl4v5.F90:1326:  write(luq_log,'(a,i4)') 'Q_ALLOCATE: klocus:',klocus
wwm_mod_xnl4v5.F90:1423:  write(luq_err,'(a,2i4)') 'Q_CHKCONFIG: iq_disp iq_geom:',iq_disp,iq_geom
wwm_mod_xnl4v5.F90:1428:  write(luq_err,'(a,i4)') 'Q_CHKCONFIG: iq_dscale:',iq_dscale
wwm_mod_xnl4v5.F90:1433:   write(luq_err,'(a,2i4)') 'Q_CHKCONFIG: iq_lump iq_gauleg:',iq_lump,iq_gauleg
wwm_mod_xnl4v5.F90:1447:  write(luq_err,'(a,2i4)') 'Q_CHKCONFIG: iq_integ:',iq_integ
wwm_mod_xnl4v5.F90:1458:  write(luq_err,'(1x,a,2i4)') 'iq_lump:',iq_lump
wwm_mod_xnl4v5.F90:1463:  write(luq_err,'(a,2i4)') 'Q_CHKCONFIG: iq_make:',iq_make
wwm_mod_xnl4v5.F90:1471:  write(luq_err,'(a)') 'Q_CHKCONFIG: Lumping or Gauss-Integration enabled when IMOD=0'
wwm_mod_xnl4v5.F90:1515:  write(luq_err,'(a,f8.2)')  'Q_CHKCONFIG: Q_DSTEP:',q_dstep
wwm_mod_xnl4v5.F90:1855:  write(luq_tst,'(a,4f11.5)') 'Q_CMPLOCUS: k1x/y k3x/y    :',k1x,k1y,k3x,k3y
wwm_mod_xnl4v5.F90:1856:  write(luq_tst,'(a,3f11.5)') 'Q_CMPLOCUS: Px Py Pmag     :',px,py,pmag
wwm_mod_xnl4v5.F90:1857:  write(luq_tst,'(a,3f11.4)') 'Q_CMPLOCUS: Pang Sang Xang :',pang*rade,sang*rade,xang*rade
wwm_mod_xnl4v5.F90:1858:  write(luq_tst,'(a,2f11.5)') 'Q_CMPLOCUS: k1m,k3m        :',k1m,k3m
wwm_mod_xnl4v5.F90:1859:  write(luq_tst,'(a,f11.5)')  'Q_CMPLOCUS: q              :',q
wwm_mod_xnl4v5.F90:1908:    write(luq_tst,'(a,4i4,2e12.4,4f12.5)') &
wwm_mod_xnl4v5.F90:1912:        write(luq_prt,'(i4,4f12.5)') iloc,x2_loc(iloc),y2_loc(iloc),&
wwm_mod_xnl4v5.F90:1920:    write(luq_err,'(a)') 'Q_CMPLOCUS: ratio > 1.5'
wwm_mod_xnl4v5.F90:1923:      write(luq_tst,'(2i5)') nlocus1,2
wwm_mod_xnl4v5.F90:1925:        write(luq_tst,'(2f12.5)') x2_loc(iloc),y2_loc(iloc)
wwm_mod_xnl4v5.F90:1941:if (iq_test >=2) write(luq_tst,'(a,4f12.5,i4)')&
wwm_mod_xnl4v5.F90:1992:    write(luq_tst,'(a,i4,4f10.5,2e12.5)') 'Q_CMPLOCUS: k2x k2y ds s jac cp:', &
wwm_mod_xnl4v5.F90:1996:if(itest >= 1) write(luq_tst,'(a,2f10.5)') 'Q_CMPLOCUS: length of locus:',sum,&
wwm_mod_xnl4v5.F90:2153:if(iq_prt>=1) write(luq_prt,'(a,i4,f16.3)') &
wwm_mod_xnl4v5.F90:2163:write(q_header,'(3i3.3,6i2.2)') naq,nkq,nlocus0,&
wwm_mod_xnl4v5.F90:2167:  write(luq_prt,'(2a)')    'Q_CTRGRID: header info:',trim(q_header)
wwm_mod_xnl4v5.F90:2168:  write(luq_prt,'(a,3i5)') 'Q_CTRGRID: naq nkq nlocus0:',naq,nkq,nlocus0
wwm_mod_xnl4v5.F90:2169:  write(luq_prt,'(a,3i3)') 'Q_CTRGRID: iq_grid,iq_geom,iq_disp:',iq_grid,iq_geom,iq_disp
wwm_mod_xnl4v5.F90:2170:  write(luq_prt,'(a,3i3)') 'Q_CTRGRID: iq_cple,iq_locus,iq_interp:',iq_cple,iq_locus,iq_interp
wwm_mod_xnl4v5.F90:2193:  if(iq_prt>=2) write(luq_prt,'(a,3f10.2,2i6)')   &
wwm_mod_xnl4v5.F90:2198:  write(bqname(5:9),'(i5.5)') min(int(q_maxdepth*10),jdep)
wwm_mod_xnl4v5.F90:2202:  write(luq_err,'(a,i4)') 'IQ_DISP=',iq_disp
wwm_mod_xnl4v5.F90:2211:if(iq_prt>=2) write(luq_prt,'(4a)')        &
wwm_mod_xnl4v5.F90:2216:  if(iq_screen>0) write(iscreen,'(2a)')   'Q_CTRGRID: Rereading of bqfile skipped: ',lastquadfile
wwm_mod_xnl4v5.F90:2217: if(iq_prt>=1) write(luq_prt,'(2a)')  'Q_CTRGRID: Rereading of bqfile skipped: ',lastquadfile
wwm_mod_xnl4v5.F90:2223:  write(luq_prt,'(2a)') 'Q_CTRGRID: Header line of grid file:',trim(q_header)
wwm_mod_xnl4v5.F90:2224:  write(luq_prt,'(2a)') 'Q_CTRGRID: Name of BINARY grid file:',trim(bqname)
wwm_mod_xnl4v5.F90:2232:if(iq_prt >= 2) write(luq_prt,'(2a,2x,2i4)') &
wwm_mod_xnl4v5.F90:2246:    write(luq_prt,'(2a)')   'Q_CTRGRID: Binary grid file detected: ',trim(bqname)
wwm_mod_xnl4v5.F90:2247:    write(luq_prt,'(a,i4)') 'Q_CTRGRID: Connected to unit:',luq_bqf
wwm_mod_xnl4v5.F90:2259:    write(luq_err,'(a)') 'BQF file deleted'
wwm_mod_xnl4v5.F90:2267:   write(luq_prt,'(a,2i4)') 'Q_CTRGRID: luq_bqf ierr:',luq_bqf,ierr
wwm_mod_xnl4v5.F90:2268:   write(luq_prt,'(2a)') 'Q_CTRGRID: r_header: ',trim(r_header)
wwm_mod_xnl4v5.F90:2272:   write(luq_prt,'(4a)') 'Q_CTRGRID: bqname  : ',trim(bqname)
wwm_mod_xnl4v5.F90:2273:   write(luq_prt,'(4a)') 'Q_CTRGRID: q_header: ',trim(q_header)
wwm_mod_xnl4v5.F90:2274:   write(luq_prt,'(4a)') 'Q_CTRGRID: r_header: ',trim(r_header)
wwm_mod_xnl4v5.F90:2285:      write(luq_prt,'(a,1x,a)') &
wwm_mod_xnl4v5.F90:2287:      write(luq_prt,'(a,1x,a)') &
wwm_mod_xnl4v5.F90:2289:      write(luq_prt,'(a)') 'Q_CTRGRID: The file headers disagree'
wwm_mod_xnl4v5.F90:2290:      write(luq_prt,'(a)') 'Q_CTRGRID: A new grid will be generated'
wwm_mod_xnl4v5.F90:2305:      write(luq_prt,'(a)')    'Q_CTRGRID: Contents of BQF file'
wwm_mod_xnl4v5.F90:2306:      write(luq_prt,'(2a)')   'Q_CTRGRID: Header:',trim(r_header)
wwm_mod_xnl4v5.F90:2307:      write(luq_prt,'(a,i4)') 'Q_CTRGRID: NK:',nkz
wwm_mod_xnl4v5.F90:2308:      write(luq_prt,'(a,i4)') 'Q_CTRGRID: NA:',naz
wwm_mod_xnl4v5.F90:2317:        write(luq_prt,'(a)') 'Q_CTRGRID: Directions do not agree'
wwm_mod_xnl4v5.F90:2319:          write(luq_prt,'(1x,a,i4,2f10.3)') 'iaz q_ad z_ad:',jaz,q_ad(jaz),z_ad(jaz)
wwm_mod_xnl4v5.F90:2331:        write(luq_prt,'(a)') 'Q_CTRGRID: Wave numbers do not agree'
wwm_mod_xnl4v5.F90:2333:          write(luq_prt,'(1x,a,i4,2f10.3)') 'ikz q_k z_sig:',jkz,q_sig(jkz),z_sig(jkz)
wwm_mod_xnl4v5.F90:2345:    write(luq_prt,'(a)') 'Q_CTRGRID: Water depths do not agree'
wwm_mod_xnl4v5.F90:2346:    write(luq_prt,'(a,2f16.2)') 'Q_CTRGRID: q_depth z_depth:',q_depth,z_depth
wwm_mod_xnl4v5.F90:2353:    if(iq_log >= 1) write(luq_log,'(a)') 'Q_CTRGRID: Existing BQF-file invalid, it will be closed'
wwm_mod_xnl4v5.F90:2377:    write(luq_log,*)
wwm_mod_xnl4v5.F90:2378:    write(luq_log,'(a)') 'Q_CTRGRID: New grid will be generated'
wwm_mod_xnl4v5.F90:2379:    write(luq_log,'(a,a)') 'Q_CTRGRID: Name of BQF file:',trim(bqname)
wwm_mod_xnl4v5.F90:2380:    write(luq_log,'(a,i4)') 'Q_CTRGRID: '//trim(bqname)//' connected to :',luq_bqf
wwm_mod_xnl4v5.F90:2383:  if(iq_screen >= 1) write(iscreen,'(2a)') &
wwm_mod_xnl4v5.F90:2400:    write(luq_log,'(a,i4)') 'Q_CTRGRID: '//trim(bqname)//' disconnected from:',luq_bqf
wwm_mod_xnl4v5.F90:2403:  if(iq_screen >=1) write(iscreen,'(a)') 'Q_CTRGRID: Grid generation completed succesfully'
wwm_mod_xnl4v5.F90:2411:  if(iq_screen >= 1) write(iscreen,'(2a)') 'Q_CTRGRID: Reading existing grid: ',trim(bqname)
wwm_mod_xnl4v5.F90:2412:  if(iq_prt >= 1)    write(luq_prt,'(2a)')  'Q_CTRGRID: Existing grid will be read:',trim(bqname)
wwm_mod_xnl4v5.F90:2413:  if(iq_log >= 1)    write(luq_log,'(2a)')  'Q_CTRGRID: Existing grid will be read:',trim(bqname)
wwm_mod_xnl4v5.F90:2436:    write(luq_err,'(a)') 'BQF file probably generated with test option off'
wwm_mod_xnl4v5.F90:2449:   write(luq_log,'(a,i4)') 'Q_CTRGRID: '//trim(bqname)//' disconnected from:',luq_bqf
wwm_mod_xnl4v5.F90:2455: if(iq_prt>=1) write(luq_prt,'(a,i4,f10.2)') &
wwm_mod_xnl4v5.F90:2460:if(iq_screen > 2) write(iscreen,'(2a,2x,f12.2)') &
wwm_mod_xnl4v5.F90:2462:if(iq_prt >=2) write(luq_prt,'(2a,2x,f12.2)') &
wwm_mod_xnl4v5.F90:2590:if(iq_test>=1) write(luq_tst,'(a,3f10.4)') 'Q_DSCALE kms,kd,q_dfac:',kms,kd,q_dfac
wwm_mod_xnl4v5.F90:2673:if(iq_log >= 1) write(luq_log,'(a,i4)') &
wwm_mod_xnl4v5.F90:2679:  write(luq_err,'(a)') q_version
wwm_mod_xnl4v5.F90:2680:  write(luq_err,'(a)')'------------------------------------------------------------'
wwm_mod_xnl4v5.F90:2687:  write(luq_err,'(a,i4)') 'Warning or non-terminating error:',iq_warn
wwm_mod_xnl4v5.F90:2688:  write(luq_err,'(a,a)')  'Name of error:',trim(err_name)
wwm_mod_xnl4v5.F90:2692:  write(luq_err,'(a,i4)') 'Terminating error:',iq_err
wwm_mod_xnl4v5.F90:2693:  write(luq_err,'(a,a)')  'Name of error:',trim(err_name)
wwm_mod_xnl4v5.F90:2694:  write(*,'(1x,a,i4)') 'Terminating error:',iq_err
wwm_mod_xnl4v5.F90:2695:  write(*,'(1x,a,a)')  'Name of error:',trim(err_name)
wwm_mod_xnl4v5.F90:2707:    if(iq_log > 0) write(luq_log,'(3a)') &
wwm_mod_xnl4v5.F90:2711:    if(iq_log >= 1) write(luq_log,'(a,i4)') &
wwm_mod_xnl4v5.F90:2724:          write(luq_err,*)
wwm_mod_xnl4v5.F90:2725:          write(luq_err,'(a)') 'Explanation of error, and recommended action'
wwm_mod_xnl4v5.F90:2726:          write(luq_err,'(a)') '--------------------------------------------'
wwm_mod_xnl4v5.F90:2727:          write(luq_err,'(a)') trim(qline)
wwm_mod_xnl4v5.F90:2739:                write(luq_err,'(a)') trim(qline)
wwm_mod_xnl4v5.F90:2754:    if(iq_log >= 1) write(luq_log,'(3a,i4)') &
wwm_mod_xnl4v5.F90:2760:  write(luq_err,*)
wwm_mod_xnl4v5.F90:2761:  write(luq_err,'(a)') 'Additional message from point of occurrence:'
wwm_mod_xnl4v5.F90:2762:  write(luq_err,'(a)') trim(err_msg)
wwm_mod_xnl4v5.F90:2763:  write(luq_err,*)
wwm_mod_xnl4v5.F90:2769:write(luq_err,'(a)') 'Trace of error'
wwm_mod_xnl4v5.F90:2770:write(luq_err,'(a)') '--------------'
wwm_mod_xnl4v5.F90:2772:  write(luq_err,'(1x,i4,2x,a)') j_stack,trim(cstack(j_stack))
wwm_mod_xnl4v5.F90:2775:write(luq_err,*)
wwm_mod_xnl4v5.F90:2930:if(iq_test >=1) write(luq_tst,'(a,6i4)') 'Q_GETLOCUS: it1,it3,itmin,iadif,ja1,ja3:',it1,it3,itmin,iadif,ja1,ja3
wwm_mod_xnl4v5.F90:2962:if(iq_test >=2) write(luq_tst,'(a,6i5)') 'Q_GETLOCUS: ik1 ia1 ik3 ia3 kmem amem:',ik1,ia1,ik3,ia3,kmem,amem
wwm_mod_xnl4v5.F90:2971:  write(luq_err,'(a,2i4)') 'Q_GETLOCUS: iamax,amem:',iamax,amem
wwm_mod_xnl4v5.F90:2982:if(iq_test>=2) write(luq_tst,'(a,i4)') 'Q_GETLOCUS: nloc:',nloc
wwm_mod_xnl4v5.F90:3158:  write(luq_trf,'(a)') '! ik1 ia1 ik3 ia3'
wwm_mod_xnl4v5.F90:3159:  write(luq_trf,'(a)') '! k1x k1y k3x k3y'
wwm_mod_xnl4v5.F90:3160:  write(luq_trf,'(a)') '! itrans kmem amem kdif iaref ibeta imirror it1 it3'
wwm_mod_xnl4v5.F90:3161:  write(luq_trf,'(a)') '! lambda depth'
wwm_mod_xnl4v5.F90:3162:  write(luq_trf,'(a)') '! nlocus'
wwm_mod_xnl4v5.F90:3163:  write(luq_trf,'(a)') '! k2x k2y k4x k4y ds jac cple sym zz'
wwm_mod_xnl4v5.F90:3165:  write(luq_trf,'(a1,2i3.3,a1,2i3.3,a1)') '(',ik1,ia1,'-',ik3,ia3,')'
wwm_mod_xnl4v5.F90:3166:  write(luq_trf,'(4f10.4)') q_k(ik1)*cos(q_a(ia1)),q_k(ik1)*sin(q_a(ia1)),&
wwm_mod_xnl4v5.F90:3168:  write(luq_trf,'(9i5)') itrans,kmem,amem,kdif,iaref,ibeta,imirror,it1,it3
wwm_mod_xnl4v5.F90:3169:  write(luq_trf,'(2f10.3)') lambda,q_depth
wwm_mod_xnl4v5.F90:3187:    write(luq_trf,'(a)') '#TRF1#'
wwm_mod_xnl4v5.F90:3188:    write(luq_trf,'(2i5)') nlocusx,8
wwm_mod_xnl4v5.F90:3190:      write(luq_trf,'(3i6,4f10.4,e13.5)') iloc,r_ik2(iloc),r_ia2(iloc),&
wwm_mod_xnl4v5.F90:3194:    write(luq_trf,'(a)') '#TRF2#'
wwm_mod_xnl4v5.F90:3195:    write(luq_trf,'(2i5)') nlocusx,8
wwm_mod_xnl4v5.F90:3197:      write(luq_trf,'(3i6,4f10.4,e13.5)') iloc,t_ik2(iloc),t_ia2(iloc),&
wwm_mod_xnl4v5.F90:3203:    write(luq_trf,'(a)') '#TRF3#'
wwm_mod_xnl4v5.F90:3204:    write(luq_trf,'(2i5)') nlocusx,9
wwm_mod_xnl4v5.F90:3219:      write(luq_trf,'(5f10.4,5e13.5)') xt2(iloc),yt2(iloc),xt4(iloc),yt4(iloc),&
wwm_mod_xnl4v5.F90:3330:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:3331:  write(luq_prt,'(a,f6.1)') 'Q_INIT:  E(f)_tail: ',qf_tail
wwm_mod_xnl4v5.F90:3332:  write(luq_prt,'(a,f6.1)') 'Q_INIT:  N(k)_tail: ',qk_tail
wwm_mod_xnl4v5.F90:3353:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:3354:  write(luq_prt,'(a)') 'Basic wave numbers, frequencies'
wwm_mod_xnl4v5.F90:3366:    write(luq_prt,'(a,i4,3f10.5,e12.4)') 'Q_INIT: ikq f sigma k k^p:', &
wwm_mod_xnl4v5.F90:3374:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:3375:  write(luq_prt,'(a)') 'Extended wave numbers and spacing'
wwm_mod_xnl4v5.F90:3405:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:3406:  write(luq_prt,'(a)') 'Q_INIT: Additional information'
wwm_mod_xnl4v5.F90:3407:  write(luq_prt,'(a,f8.1)')  'Q_depth (m):',q_depth
wwm_mod_xnl4v5.F90:3408:  write(luq_prt,'(a,i3)')    'Number of frequencies:',nkq
wwm_mod_xnl4v5.F90:3409:  write(luq_prt,'(a,f8.4)')  'Geometric f-spacing factor:',q_ffac
wwm_mod_xnl4v5.F90:3410:  write(luq_prt,'(a,f8.4)')  'Geometric k-spacing factor:',q_kfac
wwm_mod_xnl4v5.F90:3411:  write(luq_prt,'(a,2f8.3)') 'fmin fmax (Hz):',fqmin,fqmax
wwm_mod_xnl4v5.F90:3412:  write(luq_prt,'(a,2f8.3)') 'kmin kmax (Hz):',kqmin,kqmax
wwm_mod_xnl4v5.F90:3413:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:3415:  write(luq_prt,*) '     i      f         df       sig      dsig       k         dk         cg'
wwm_mod_xnl4v5.F90:3418:    write(luq_prt,'(1x,i4,7f10.4)') &
wwm_mod_xnl4v5.F90:3436:if(iq_prt >= 2) write(luq_prt,'(a,i4)') &
wwm_mod_xnl4v5.F90:3452:if(iq_prt >= 2) write(luq_prt,'(a,2i4)') &
wwm_mod_xnl4v5.F90:3463:  if(iq_prt>0) write(luq_prt,'(a)') 'Q_INIT: take care of q_dird1 and check if sector is OK'
wwm_mod_xnl4v5.F90:3475:  write(luq_prt,'(a,3f10.3)')     'Q_INIT: d(1),d(n),dsector:',q_dird1,q_dird2,q_sector
wwm_mod_xnl4v5.F90:3476:  write(luq_prt,'(a,f6.2,a)')     'Q_INIT: Angular step     :',q_deltad,' degrees'
wwm_mod_xnl4v5.F90:3477:  write(luq_prt,'(a,2f8.2,i4,a)') 'Q_INIT: ang1 ang2 nang   :',q_ang1,q_ang2,naq,' degrees'
wwm_mod_xnl4v5.F90:3478:  write(luq_prt,'(a,i4)')         'Q_INIT: #Angles on circle:',ncirc
wwm_mod_xnl4v5.F90:3479:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:3488:    write(luq_prt,'(a,i4,f10.4,f10.2)') 'Q_INIT: iaq q_a q_ad:',iaq,q_a(iaq),q_ad(iaq)
wwm_mod_xnl4v5.F90:3489:    if(iaq==naq) write(luq_prt,*)
wwm_mod_xnl4v5.F90:3510:  write(luq_tst,'(a,3i4)') 'Q_INIT: iq_grid iaref iamax:',iq_grid,iaref,iamax
wwm_mod_xnl4v5.F90:3511:  write(luq_tst,'(a,4i4)') 'Q_INIT: iaq1 iaq2 iag1 iag2:',iaq1,iaq2,iag1,iag2
wwm_mod_xnl4v5.F90:3515:  write(luq_trf,'(a)') '#GRIDINFO#'
wwm_mod_xnl4v5.F90:3516:  write(luq_trf,'(2i4)') nkq,naq
wwm_mod_xnl4v5.F90:3517:  write(luq_trf,'(10f8.4)') q_k
wwm_mod_xnl4v5.F90:3518:  write(luq_trf,'(10f8.2)') q_a*rade
wwm_mod_xnl4v5.F90:3672:  if(itest >= 1) write(luq_tst,'(a,6e12.5)') &
wwm_mod_xnl4v5.F90:3694:     if(itest>=2) write(luq_tst,'(a,i4)') 'Q_LOCPOS/Z_ROOT2/IERR/1=',ierr
wwm_mod_xnl4v5.F90:3697:    if(itest >= 1) write(luq_tst,'(a,4f12.5)') &
wwm_mod_xnl4v5.F90:3714:      if(itest >= 2) write(luq_tst,'(a,i4,3f12.5,2e13.5)') &
wwm_mod_xnl4v5.F90:3727:   if(ierr>0 .and. itest>=2) write(luq_tst,'(a,i4)') 'Q_LOCPOS/Z_ROOT2/IERR/2=',ierr
wwm_mod_xnl4v5.F90:3748:     if(itest >= 2) write(luq_tst,'(a,i4,3f12.5,2e12.5)') &
wwm_mod_xnl4v5.F90:3761:   if(ierr > 0 .and. itest>=2) write(luq_tst,'(a,i4)') 'Q_LOCPOS/Z_ROOT2/IERR/3=',ierr
wwm_mod_xnl4v5.F90:3777:      if(itest >=  2) write(luq_tst,'(a,i4,3f12.5,2e12.5)') &
wwm_mod_xnl4v5.F90:3790:   if(ierr>0.and.itest>=2) write(luq_tst,'(a,i4)') 'Q_LOCPOS/Z_ROOT2/IERR/4=',ierr
wwm_mod_xnl4v5.F90:3799:    if(itest >= 1) write(luq_tst,'(a,6e12.5)') &
wwm_mod_xnl4v5.F90:3816:if(itest >= 1) write(luq_tst,'(a,3f12.6)') &
wwm_mod_xnl4v5.F90:3832:if(itest >= 1) write(luq_tst,'(a,4f10.5,2e12.5)') 'Q_LOCPOS: k1 k2 z1/2:',&
wwm_mod_xnl4v5.F90:3842:  if(itest >= 1) write(luq_tst,'(a,i4,4f12.5,2e12.5)') &
wwm_mod_xnl4v5.F90:3850:  write(luq_tst,'(a,2f10.4,2e13.5)') 'Q_LOCPOS: beta1/2 xlocus2(beta1/2):',beta1,beta2,x_locus2(beta1),x_locus2(beta2)
wwm_mod_xnl4v5.F90:3857:  if(itest>=2) write(luq_tst,'(a,i4)') 'Q_LOCPOS/Z_ROOT2/IERR_W/5=',ierr
wwm_mod_xnl4v5.F90:3867:if(itest >= 1) write(luq_tst,'(a,4f12.6,e12.5)') &
wwm_mod_xnl4v5.F90:3893:  write(luq_tst,'(a,4f10.5)') 'Q_LOCPOS: aa,bb,mm,mm1:',aa,bb,mm,mm1
wwm_mod_xnl4v5.F90:3894:  write(luq_tst,'(a,4f10.5)') 'Q_LOCPOS: length of ellipse:',loclen
wwm_mod_xnl4v5.F90:4045:  if(iq_screen==2) write(iscreen,*) 'k1-ring:',ikq1
wwm_mod_xnl4v5.F90:4055: if(iq_test >=1) write(luq_tst,'(a,i4,2x,2f8.4)') &
wwm_mod_xnl4v5.F90:4059:   if(iq_screen==2) write(iscreen,*) 'k1-k3 indices:',ikq1,ikq3
wwm_mod_xnl4v5.F90:4063:   if(iq_test >= 1) write(luq_tst,'(a,3f12.6)') &
wwm_mod_xnl4v5.F90:4077:        write(luq_tst,'(a,2f10.4,2i4)') 'Q_MAKEGRID: k3 a3 ikq3 iaq3: ',kk3,aa3,ikq3,iaq3
wwm_mod_xnl4v5.F90:4078:        write(luq_tst,'(a,4f11.5)')     'Q_MAKEGRID: k1x/y k3x/y    :',k1x,k1y,k3x,k3y
wwm_mod_xnl4v5.F90:4222:!     write(luq_prt,'(a,4i5)') 'Q_MAKEGRID kmem amem nlocus:',kmem,amem,nlocus,nzloc
wwm_mod_xnl4v5.F90:4231:write(luq_bqf) q_header
wwm_mod_xnl4v5.F90:4237:write(luq_bqf) naq,nkq
wwm_mod_xnl4v5.F90:4238:write(luq_bqf) q_sig
wwm_mod_xnl4v5.F90:4239:write(luq_bqf) q_ad
wwm_mod_xnl4v5.F90:4240:write(luq_bqf) iq_geom,iq_disp,iq_geom
wwm_mod_xnl4v5.F90:4241:write(luq_bqf) q_depth
wwm_mod_xnl4v5.F90:4247:write(luq_bqf) quad_nloc
wwm_mod_xnl4v5.F90:4248:write(luq_bqf) quad_ik2
wwm_mod_xnl4v5.F90:4249:write(luq_bqf) quad_ia2
wwm_mod_xnl4v5.F90:4250:write(luq_bqf) quad_ik4
wwm_mod_xnl4v5.F90:4251:write(luq_bqf) quad_ia4
wwm_mod_xnl4v5.F90:4252:write(luq_bqf) quad_w1k2
wwm_mod_xnl4v5.F90:4253:write(luq_bqf) quad_w2k2
wwm_mod_xnl4v5.F90:4254:write(luq_bqf) quad_w3k2
wwm_mod_xnl4v5.F90:4255:write(luq_bqf) quad_w4k2
wwm_mod_xnl4v5.F90:4256:write(luq_bqf) quad_w1k4
wwm_mod_xnl4v5.F90:4257:write(luq_bqf) quad_w2k4
wwm_mod_xnl4v5.F90:4258:write(luq_bqf) quad_w3k4
wwm_mod_xnl4v5.F90:4259:write(luq_bqf) quad_w4k4
wwm_mod_xnl4v5.F90:4260:write(luq_bqf) quad_zz
wwm_mod_xnl4v5.F90:4261:write(luq_bqf) quad_t2
wwm_mod_xnl4v5.F90:4262:write(luq_bqf) quad_t4
wwm_mod_xnl4v5.F90:4264:write(luq_bqf) quad_jac
wwm_mod_xnl4v5.F90:4265:write(luq_bqf) quad_cple
wwm_mod_xnl4v5.F90:4266:write(luq_bqf) quad_sym
wwm_mod_xnl4v5.F90:4267:write(luq_bqf) quad_ws
wwm_mod_xnl4v5.F90:4272:if(iq_screen >= 1 .and. iq_test>=1) write(iscreen,'(2a)') 'Q_MAKEGRID: LASTQUADFILE: ',lastquadfile
wwm_mod_xnl4v5.F90:4284:    write(luq_log,*)
wwm_mod_xnl4v5.F90:4285:    write(luq_log,'(5a)') 'Q_MAKEGRID: Grid files ',trim(aqname),' and ',trim(bqname),' deleted'
wwm_mod_xnl4v5.F90:4286:    write(luq_log,'(a)') 'Q_MAKEGRID: Since an error occurred during the generation'
wwm_mod_xnl4v5.F90:4287:    write(luq_log,'(a)') 'Q_MAKEGRID: of the interaction grid'
wwm_mod_xnl4v5.F90:4295:  write(luq_prt,'(a,i10)') 'Q_MAKEGRID: Total number of points on loci        :',nztot2
wwm_mod_xnl4v5.F90:4296:  write(luq_prt,'(a,i10)') 'Q_MAKEGRID: Total number of stored points on locus:',nztot1
wwm_mod_xnl4v5.F90:4297:  write(luq_prt,'(a,i10)') 'Q_MAKEGRID: Total number of zero points on locus  :',nztot2-nztot1
wwm_mod_xnl4v5.F90:4298:  write(luq_prt,'(a,f8.2)') 'Q_MAKEGRID: Reduction factor (%):',real(nztot2-nztot1)/real(nztot2)*100.
wwm_mod_xnl4v5.F90:4413:  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_mod   :',iq_mod
wwm_mod_xnl4v5.F90:4414:  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_xdia  :',iq_xdia
wwm_mod_xnl4v5.F90:4415:  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_lump  :',iq_lump
wwm_mod_xnl4v5.F90:4416:  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_gauleg:',iq_gauleg
wwm_mod_xnl4v5.F90:4456:  if(itest>=1) write(luq_tst,'(a,2i4)') 'Q_MODIFY nold nnew:',nlocus1,nnew
wwm_mod_xnl4v5.F90:4497:      write(luq_tst,'(a,2f10.4,i4)') 'Q_MODIFY: GAULEG x1,x2,n:',zero,slen,nnew
wwm_mod_xnl4v5.F90:4498:      write(luq_tst,'(a)') 'Q_MODIFY: Gauss-Legendre spacing'
wwm_mod_xnl4v5.F90:4499:      write(luq_tst,'(10f12.4)') (snew(inew),inew=1,nnew)
wwm_mod_xnl4v5.F90:4500:      write(luq_tst,'(a)') 'Q_MODIFY: Gauss-Legendre weights'
wwm_mod_xnl4v5.F90:4501:      write(luq_tst,'(10f12.4)') (ds_mod(inew),inew=1,nnew)
wwm_mod_xnl4v5.F90:4519:    write(luq_tst,'(a,2f12.5)') 'Q_MODIFY: Slen q:',slen,q
wwm_mod_xnl4v5.F90:4520:    write(luq_tst,'(a,i4)') 'Q_MODIFY: nold /sold:',nold
wwm_mod_xnl4v5.F90:4521:    write(luq_tst,'(10f12.6)') sold
wwm_mod_xnl4v5.F90:4522:    write(luq_tst,'(a,i4)') 'Q_MODIFY: nnew /snew:',nnew
wwm_mod_xnl4v5.F90:4523:    write(luq_tst,'(10f12.6)') snew
wwm_mod_xnl4v5.F90:4524:    write(luq_tst,'(a)') 'Q_MODIFY: x2_loc'
wwm_mod_xnl4v5.F90:4525:    write(luq_tst,'(10f13.5)') (x2_loc(iloc), iloc=1,nold)
wwm_mod_xnl4v5.F90:4526:    write(luq_tst,'(a)') 'Q_MODIFY: y2_loc'
wwm_mod_xnl4v5.F90:4527:    write(luq_tst,'(10f13.5)') (y2_loc(iloc), iloc=1,nold)
wwm_mod_xnl4v5.F90:4537:    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4541:    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4545:    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4549:    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4553:    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4569:        if(itest>=1) write(luq_tst,'(a,2i4,f8.3,3e12.4,f4.0,e12.4)') &
wwm_mod_xnl4v5.F90:4582:      if(ierr > 0) write(luq_err,*) 'Z_INTP1 jac_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4586:      if(ierr > 0) write(luq_err,*) 'Z_INTP1 cp_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4594:    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4598:    if(ierr > 0) write(luq_err,*) 'Z_INTP1 y_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4602:    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4606:    if(ierr > 0) write(luq_err,*) 'Z_INTP1 y_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4610:    if(ierr > 0) write(luq_err,*) 'Z_INTP1 s_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4614:    write(luq_tst,'(a)') 'Q_MODIFY: s_loc'
wwm_mod_xnl4v5.F90:4615:    write(luq_tst,'(10f13.5)') (s_loc(iloc), iloc=1,nold)
wwm_mod_xnl4v5.F90:4616:    write(luq_tst,'(a)') 'Q_MODIFY: s_mod'
wwm_mod_xnl4v5.F90:4617:    write(luq_tst,'(10f13.5)') (s_mod(iloc), iloc=1,nold)
wwm_mod_xnl4v5.F90:4634:        if(itest>=1) write(luq_tst,'(a,2i4,f8.3,3e12.4,f4.0,e12.4)') &
wwm_mod_xnl4v5.F90:4647:      if(ierr > 0) write(luq_err,*) 'Z_INTP1 jac_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4651:      if(ierr > 0) write(luq_err,*) 'Z_INTP1 cp_loc, ierr=',ierr
wwm_mod_xnl4v5.F90:4669:  write(luq_tst,'(a)') 'Q_MODIFY: x2_mod'
wwm_mod_xnl4v5.F90:4670:  write(luq_tst,'(10f12.5)') (x2_mod(iloc),iloc=1,nlocus)
wwm_mod_xnl4v5.F90:4671:  write(luq_tst,'(a)') 'Q_MODIFY: y2_mod'
wwm_mod_xnl4v5.F90:4672:  write(luq_tst,'(10f12.5)') (y2_mod(iloc),iloc=1,nlocus)
wwm_mod_xnl4v5.F90:4673:  write(luq_tst,'(a)') 'Q_MODIFY: x4_mod'
wwm_mod_xnl4v5.F90:4674:  write(luq_tst,'(10f12.5)') (x4_mod(iloc),iloc=1,nlocus)
wwm_mod_xnl4v5.F90:4675:  write(luq_tst,'(a)') 'Q_MODIFY: y4_mod'
wwm_mod_xnl4v5.F90:4676:  write(luq_tst,'(10f12.5)') (y4_mod(iloc),iloc=1,nlocus)
wwm_mod_xnl4v5.F90:4677:  write(luq_tst,'(a)') 'Q_MODIFY: s_mod'
wwm_mod_xnl4v5.F90:4678:  write(luq_tst,'(10f12.5)') (s_mod(iloc),iloc=1,nlocus)
wwm_mod_xnl4v5.F90:4679:  write(luq_tst,'(a)') 'Q_MODIFY: ds_loc'
wwm_mod_xnl4v5.F90:4680:  write(luq_tst,'(10f12.5)') (ds_loc(iloc),iloc=1,nold)
wwm_mod_xnl4v5.F90:4681:  write(luq_tst,'(a)') 'Q_MODIFY: ds_mod'
wwm_mod_xnl4v5.F90:4682:  write(luq_tst,'(10f12.5)') (ds_mod(iloc),iloc=1,nlocus)
wwm_mod_xnl4v5.F90:4712:  write(luq_tst,'(a)') 'Q_MODIFY: k2m_mod'
wwm_mod_xnl4v5.F90:4713:  write(luq_tst,'(10f12.5)') (k2m_mod(iloc),iloc=1,nlocus)
wwm_mod_xnl4v5.F90:4714:  write(luq_tst,'(a)') 'Q_MODIFY: k2a_mod'
wwm_mod_xnl4v5.F90:4715:  write(luq_tst,'(10f12.5)') (k2a_mod(iloc),iloc=1,nlocus)
wwm_mod_xnl4v5.F90:4716:  write(luq_tst,'(a)') 'Q_MODIFY: k4m_mod'
wwm_mod_xnl4v5.F90:4717:  write(luq_tst,'(10f12.5)') (k4m_mod(iloc),iloc=1,nlocus)
wwm_mod_xnl4v5.F90:4718:  write(luq_tst,'(a)') 'Q_MODIFY: k4a_mod'
wwm_mod_xnl4v5.F90:4719:  write(luq_tst,'(10f12.5)') (k4a_mod(iloc),iloc=1,nlocus)
wwm_mod_xnl4v5.F90:4720:  write(luq_tst,'(a)') 'Q_MODIFY: sym_mod'
wwm_mod_xnl4v5.F90:4721:  write(luq_tst,'(20f3.0)') (sym_mod(iloc),iloc=1,nlocus)
wwm_mod_xnl4v5.F90:4809:if(iq_test>=2) write(luq_tst,'(a,2i4,4f10.5)') 'Q_NEAREST-A:',ik,ia,w1,w2,w3,w4
wwm_mod_xnl4v5.F90:4863:  write(luq_tst,'(a,2i4,4f10.5)') 'Q_NEAREST-B:',ik,ia,w1,w2,w3,w4
wwm_mod_xnl4v5.F90:4864:  write(luq_tst,*)
wwm_mod_xnl4v5.F90:5031:  if(iq_test>=2) write(luq_tst,'(a,3f12.6)') 'Q_POLAR2: loclen dsz dk:',loclen,dsz,dk
wwm_mod_xnl4v5.F90:5075:  if(iq_test>=2) write(luq_tst,'(a,i4,3f11.6)') 'Q_POLAR2: npol kmin kmax kratio:',npol,kmin,kmax,kratio
wwm_mod_xnl4v5.F90:5110:  write(luq_tst,'(a,3i4)')           'Q_POLAR2: nlocus0 npol nlocus1:',nlocus0,npol,nlocus1
wwm_mod_xnl4v5.F90:5111:  write(luq_tst,'(a,2f12.6,i4)')     'Q_POLAR2: kmin kmax iq_locus  :',kmin,kmax,iq_locus
wwm_mod_xnl4v5.F90:5112:  if(iq_locus==1) write(luq_tst,'(a,f10.4)') 'Q_POLAR2: dk                  :',dk
wwm_mod_xnl4v5.F90:5113:  if(iq_locus==3) write(luq_tst,'(a,f10.4)') 'Q_POLAR2: kratio              :',kratio
wwm_mod_xnl4v5.F90:5115:    write(luq_tst,'(a,i4,4f13.7)') &
wwm_mod_xnl4v5.F90:5284:  if(iq_screen>0) write(iscreen,'(a,i4)') 'Q_SETCONFIG: iquad=',iquad
wwm_mod_xnl4v5.F90:5286:  write(luq_err,'(a,i4)') 'Q_SETCONFIG: Value of IQUAD:',iquad
wwm_mod_xnl4v5.F90:5311:    write(luq_log,*)
wwm_mod_xnl4v5.F90:5312:    write(luq_log,'(a)') 'Q_SETCONFIG: Configuration file '//trim(qbase)//'.cfg has been found'
wwm_mod_xnl4v5.F90:5313:    write(luq_log,'(a,i4)') 'Q_SETCONFIG: '//trim(qbase)//'.cfg connected to :',luq_cfg
wwm_mod_xnl4v5.F90:5340:          if(iq_screen>0) write(iscreen,'(a)') 'Q_SETCONFIG: geometric scaling disabled'
wwm_mod_xnl4v5.F90:5341:          if(iq_prt>=1)   write(luq_prt,'(a)') 'Q_SETCONFIG: geometric scaling disabled'
wwm_mod_xnl4v5.F90:5373:  if(iq_log >= 1) write(luq_log,'(a,i4)') &
wwm_mod_xnl4v5.F90:5379:    write(luq_log,*)
wwm_mod_xnl4v5.F90:5380:    write(luq_log,'(a)') 'Q_SETCONFIG: Configuration file '//trim(qbase)//'.CFG has not been found'
wwm_mod_xnl4v5.F90:5483:if(iq_prt>=1) write(luq_prt,'(a,f10.2)') 'Q_SEARCHGRID: Input target depth:',depth
wwm_mod_xnl4v5.F90:5489:if(iq_prt>=1) write(luq_prt,'(a,i4,f12.2)') &
wwm_mod_xnl4v5.F90:5494:     write(luq_prt,'(a,f10.2)') 'Q_SEARCHGRID: target depth:',q_depth
wwm_mod_xnl4v5.F90:5495:     write(iscreen,'(a)') 'Q_SEARCHGRID: grid accepted, read whole database'
wwm_mod_xnl4v5.F90:5497:  if(iq_screen>=1) write(iscreen,'(a)') 'Q_SEARCHGRID: grid accepted, read whole database'
wwm_mod_xnl4v5.F90:5514:if(iq_prt>=2) write(luq_prt,'(a,3i6)') 'Q_SEARCHGRID: idepth,id_lower/upper:',&
wwm_mod_xnl4v5.F90:5522:  if(iq_prt>=2) write(luq_prt,'(a,i6,f8.1)') 'Q_SEARCHGRID: downwards  id q_depth:',id,q_depth
wwm_mod_xnl4v5.F90:5526:  if(iq_prt>=2) write(luq_prt,'(a,i4)') 'Q_SEARCHGRID: igrid:',igrid
wwm_mod_xnl4v5.F90:5529:    if(iq_prt>=2) write(luq_prt,'(a,f8.2)') 'Q_SEARCHGRID: valid grid found for depth:',q_depth
wwm_mod_xnl4v5.F90:5542:  if(iq_prt>=2) write(luq_prt,'(a,i6,f8.1)') 'Q_SEARCHGRID: upwards  id q_depth:',id,q_depth
wwm_mod_xnl4v5.F90:5546:  if(iq_prt>=2) write(luq_prt,'(a,i4)') 'Q_SEARCHGRID: igrid:',igrid
wwm_mod_xnl4v5.F90:5549:    if(iq_prt>=2) write(luq_prt,'(a,f8.2)') 'Q_SEARCHGRID: valid grid found for depth:',q_depth
wwm_mod_xnl4v5.F90:5554:if(iq_prt>=1) write(luq_prt,*)
wwm_mod_xnl4v5.F90:5573:  write(luq_prt,'(a,3f8.2)') 'Q_SEARCHGRID: d_lower d_target d_upper      :',d_lower,s_depth,d_upper
wwm_mod_xnl4v5.F90:5574:  write(luq_prt,'(a,2f8.2)') 'Q_SEARCHGRID: r_lower r_upper               :',r_lower,r_upper
wwm_mod_xnl4v5.F90:5595:if(iq_prt>=1) write(luq_prt,'(a,2f10.2)') &
wwm_mod_xnl4v5.F90:5597:if(iq_screen>0) write(iscreen,'(a,f10.2)') &
wwm_mod_xnl4v5.F90:5609:  write(luq_prt,'(a,2f8.4)') 'Q_SEARCHGRID: target and nearest scale factors:',dfac1,dfac2
wwm_mod_xnl4v5.F90:5610:  write(luq_prt,'(a,f8.4)')  'Q_SEARCHGRID: compound scale factor           :',q_scale
wwm_mod_xnl4v5.F90:5617:  write(luq_prt,'(a,f12.2)') 'Q_SEARCHGRID: Q_CTRGRID called with depth:',q_depth
wwm_mod_xnl4v5.F90:5618:  write(luq_prt,'(a,i4)') 'Q_SEARCHGRID: igrid of nearest grid operation:',igrid
wwm_mod_xnl4v5.F90:5626:write(*,*) 'q_searchgrid q_depth on exit:',q_depth
wwm_mod_xnl4v5.F90:5712:  if(iq_prt>0)  write(luq_prt,'(2a)') 'TRACE -> ',trim(mod_name)
wwm_mod_xnl4v5.F90:5713:  if(iq_test>0) write(luq_tst,'(2a)') 'TRACE -> ',trim(mod_name)
wwm_mod_xnl4v5.F90:5714:  if(iq_screen >= 2) write(iscreen,'(2a)') 'TRACE -> ',trim(mod_name)
wwm_mod_xnl4v5.F90:5742:    write(luq_err,'(a)') 'Module name:',mod_name
wwm_mod_xnl4v5.F90:5827:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5828:  write(luq_prt,'(a)') 'Summary of settings for QUAD computation'
wwm_mod_xnl4v5.F90:5829:  write(luq_prt,'(a)') '------------------------------------------------'
wwm_mod_xnl4v5.F90:5830:  write(luq_prt,'(a,i4)')     'Number of wave numbers          :',nkq
wwm_mod_xnl4v5.F90:5831:  write(luq_prt,'(a,i4)')     'Number of directions            :',naq
wwm_mod_xnl4v5.F90:5832:  write(luq_prt,'(a,f10.5)')  'Minimum frequency (Hz)          :',fqmin
wwm_mod_xnl4v5.F90:5833:  write(luq_prt,'(a,f10.5)')  'Maximum frequency (Hz)          :',fqmax
wwm_mod_xnl4v5.F90:5834:  write(luq_prt,'(a,f10.2)')  'Water depth (m)                 :',q_depth
wwm_mod_xnl4v5.F90:5835:  write(luq_prt,'(a,i4)')     'Preferred number of locus points:',nlocus0
wwm_mod_xnl4v5.F90:5837:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5838:  write(luq_prt,'(a,f10.3)') 'Gravitational acceleration:',q_grav
wwm_mod_xnl4v5.F90:5839:!  write(luq_prt,'(a,f10.3)') 'Density of water          :',q_rhow
wwm_mod_xnl4v5.F90:5840:!  write(luq_prt,'(a,f10.2)') 'Power spectral tail E(f)  :',qf_tail
wwm_mod_xnl4v5.F90:5841:!  write(luq_prt,'(a,f10.2)') 'Power spectral tail N(k)  :',qk_tail
wwm_mod_xnl4v5.F90:5843:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5844:  if(iq_type==1) write(luq_prt,'(a)') 'IQUAD = 1: Deep water'
wwm_mod_xnl4v5.F90:5845:  if(iq_type==2) write(luq_prt,'(a)') 'IQUAD = 2: Deep water & WAM depth scaling'
wwm_mod_xnl4v5.F90:5846:  if(iq_type==3) write(luq_prt,'(a)') 'IQUAD = 3: Direct finite depth calculation'
wwm_mod_xnl4v5.F90:5847:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5849:  write(luq_prt,'(a,f5.2)') 'Step size in m of BQF coding:',q_dstep
wwm_mod_xnl4v5.F90:5850:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5852:  if(iq_grid==1) write(luq_prt,'(a)') 'Symmetric sector grid'
wwm_mod_xnl4v5.F90:5853:  if(iq_grid==2) write(luq_prt,'(a)') 'Non-symmetric sector grid'
wwm_mod_xnl4v5.F90:5854:  if(iq_grid==3) write(luq_prt,'(a)') 'Non-symmetric full circle grid'
wwm_mod_xnl4v5.F90:5856:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5857:  if(iq_compact==0) write(luq_prt,'(a)') 'No compacting of data along locus'
wwm_mod_xnl4v5.F90:5858:  if(iq_compact==1) write(luq_prt,'(a)') 'Compact data along locus by eliminating zero contributions'
wwm_mod_xnl4v5.F90:5860:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5861:  if(iq_dscale==0) write(luq_prt,'(a)') 'No WAM depth scaling'
wwm_mod_xnl4v5.F90:5862:  if(iq_dscale==1) write(luq_prt,'(a)') 'WAM depth scaling of transfer'
wwm_mod_xnl4v5.F90:5864:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5865:  if(iq_screen==0) write(luq_prt,'(a)') 'No output to screen'
wwm_mod_xnl4v5.F90:5866:  if(iq_screen>=1) write(luq_prt,'(a)') 'Intermediate output to screen'
wwm_mod_xnl4v5.F90:5867:  if(iq_screen>=2) write(luq_prt,'(a)') 'Intermediate output to screen + subroutine tracing'
wwm_mod_xnl4v5.F90:5868:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5870:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5871:  if(iq_search==0) write(luq_prt,'(a)') 'No search is carried out for nearest QUAD grid'
wwm_mod_xnl4v5.F90:5872:  if(iq_search==1) write(luq_prt,'(a)') 'A search is carried out for nearest QUAD grid'
wwm_mod_xnl4v5.F90:5874:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5875:  if(iq_gauleg==0) write(luq_prt,'(a)')    'Rectangular integration'
wwm_mod_xnl4v5.F90:5876:  if(iq_gauleg>0)  write(luq_prt,'(a,i4)') 'Gauss-Legendre integration with N=',iq_gauleg
wwm_mod_xnl4v5.F90:5878:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5879:  if(iq_cple==1) write(luq_prt,'(a)') 'Deep water coupling coefficient of Webb'
wwm_mod_xnl4v5.F90:5880:  if(iq_cple==2) write(luq_prt,'(a)') 'Finite depth coupling coefficient of H&H'
wwm_mod_xnl4v5.F90:5881:  if(iq_cple==3) write(luq_prt,'(a)') 'Finite depth coupling coefficient of Gorman'
wwm_mod_xnl4v5.F90:5882:  if(iq_cple==4) write(luq_prt,'(a)') 'Deep water coefficient of Zakharov'
wwm_mod_xnl4v5.F90:5883:  if(iq_cple==5) write(luq_prt,'(a)') 'Finite depth coefficient of Zakharov'
wwm_mod_xnl4v5.F90:5885:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5886:  if(iq_disp==1) write(luq_prt,'(a)') 'Deep water dispersion relation'
wwm_mod_xnl4v5.F90:5887:  if(iq_disp==2) write(luq_prt,'(a)') 'Finite depth linear dispersion relation'
wwm_mod_xnl4v5.F90:5888:  if(iq_disp==3) write(luq_prt,'(a)') 'Non linear finite depth dispersion'
wwm_mod_xnl4v5.F90:5890:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5891:  if(iq_filt==0) write(luq_prt,'(a)') 'Filtering of quadruplets off'
wwm_mod_xnl4v5.F90:5893:     write(luq_prt,'(a)') 'Filtering of quadruplets on'
wwm_mod_xnl4v5.F90:5894:     write(luq_prt,*)
wwm_mod_xnl4v5.F90:5895:     write(luq_prt,'(a,f8.2)')  'Maximum ratio of k1 and k3        :',qf_krat
wwm_mod_xnl4v5.F90:5896:     write(luq_prt,'(a,f8.2)')  'Maximum directional difference    :',qf_dmax
wwm_mod_xnl4v5.F90:5897:     write(luq_prt,'(a,e12.3)') 'Fraction of maximum energy density:',qf_frac
wwm_mod_xnl4v5.F90:5900:!  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5901:!  if(iq_geom==0) write(luq_prt,'(a)') 'Only directional scaling of loci'
wwm_mod_xnl4v5.F90:5902:!  if(iq_geom==1) write(luq_prt,'(a)') 'Geometric scaling of loci using R-T method'
wwm_mod_xnl4v5.F90:5904:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5905:  if(iq_locus==1) write(luq_prt,'(a)') 'Compute locus with polar method with fixed k-step'
wwm_mod_xnl4v5.F90:5906:  if(iq_locus==2) write(luq_prt,'(a)') 'Compute locus with polar method using adaptive k-step'
wwm_mod_xnl4v5.F90:5907:  if(iq_locus==3) write(luq_prt,'(a)') 'Compute locus with polar method using geometric k-step'
wwm_mod_xnl4v5.F90:5909:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5910:  if(iq_sym==0)  write(luq_prt,'(a)') 'Handling of symmetries disabled'
wwm_mod_xnl4v5.F90:5911:  if(iq_sym==1)  write(luq_prt,'(a)') 'Handling of symmetries enabled'
wwm_mod_xnl4v5.F90:5913:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5914:  if(iq_make==1) write(luq_prt,'(a)') 'Make quadruplet grid when necessary'
wwm_mod_xnl4v5.F90:5915:  if(iq_make==2) write(luq_prt,'(a)') 'Always make quadruplet grid'
wwm_mod_xnl4v5.F90:5916:  if(iq_make==3) write(luq_prt,'(a)') 'Stop after generation of quadruplet grid'
wwm_mod_xnl4v5.F90:5918:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5919:  if(iq_interp==1) write(luq_prt,'(a)') 'Apply bi-linear interpotion to retrieve action density'
wwm_mod_xnl4v5.F90:5920:  if(iq_interp==2) write(luq_prt,'(a)') 'Take nearest bin to retrieve action density'
wwm_mod_xnl4v5.F90:5922:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5923:  if(iq_lump==0) write(luq_prt,'(a)') 'Lumping of coefficients along locus disabled'
wwm_mod_xnl4v5.F90:5924:  if(iq_lump>0)  write(luq_prt,'(a)') 'Lumping of coefficients along locus enabled'
wwm_mod_xnl4v5.F90:5926:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5927:  if(iq_mod==0) write(luq_prt,'(a)') '?X? Spacing of point along locus as initially computed'
wwm_mod_xnl4v5.F90:5928:  if(iq_mod==1) write(luq_prt,'(a)') 'Equidistant spacing of points along locus'
wwm_mod_xnl4v5.F90:5930:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5931:  if(iq_tail==0) write(luq_prt,'(a)') 'No parametric tail is added'
wwm_mod_xnl4v5.F90:5932:  if(iq_tail==1) write(luq_prt,'(a,f8.2,a)') 'Parametric tail is added from ',ff_tail,' times maximum frequency'
wwm_mod_xnl4v5.F90:5934:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5935:  if(iq_trace==0) write(luq_prt,'(a)') 'Subroutine tracing disabled'
wwm_mod_xnl4v5.F90:5936:  if(iq_trace>0)  write(luq_prt,'(a)') 'Subroutine tracing enabled'
wwm_mod_xnl4v5.F90:5939:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5940:  write(luq_prt,'(a,i4)') 'IQ_INTEG:',iq_integ
wwm_mod_xnl4v5.F90:5941:  if(iq_integ==0) write(luq_prt,'(a)') 'No test output of integration'
wwm_mod_xnl4v5.F90:5942:  if(iq_integ==1) write(luq_prt,'(a)') 'Summary output of integration per locus'
wwm_mod_xnl4v5.F90:5943:  if(iq_integ==2) write(luq_prt,'(a)') 'Extended output of integration along locus'
wwm_mod_xnl4v5.F90:5944:  if(iq_integ==3) write(luq_prt,'(a)') 'Line function along locus'
wwm_mod_xnl4v5.F90:5946:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5947:  if(iq_t13==0) write(luq_prt,'(a)') 'No test output of T13 integration'
wwm_mod_xnl4v5.F90:5948:  if(iq_t13==1) write(luq_prt,'(a)') 'Summary output of T13 integration'
wwm_mod_xnl4v5.F90:5950:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5953:!    write(luqprt,'(a)') 'Start point for locus according to Resio&Tracy'
wwm_mod_xnl4v5.F90:5955:!    write(luqprt,'(a)') 'Start point for locus equal to k3'
wwm_mod_xnl4v5.F90:5957:  write(luq_prt,*)
wwm_mod_xnl4v5.F90:5958:  write(luq_prt,'(a,i4)') 'Level of printed output        :',iq_prt
wwm_mod_xnl4v5.F90:5959:  write(luq_prt,'(a,i4)') 'Level of logging output        :',iq_log
wwm_mod_xnl4v5.F90:5960:  write(luq_prt,'(a,i4)') 'Level of test output           :',iq_test
wwm_mod_xnl4v5.F90:5961:  write(luq_prt,'(a,i4)') 'Level of trace output          :',iq_trace
wwm_mod_xnl4v5.F90:5962:  write(luq_prt,'(a,i4)') 'Level of transformation output :',iq_trf
wwm_mod_xnl4v5.F90:5963:  write(luq_prt,'(a)')   '----------------------------------------------'
wwm_mod_xnl4v5.F90:6177:if(iq_test >=2) write(luq_tst,'(a,4i3)') 'Q_T13V4: ik1,ia1 ik3 ia3:',ik1,ia1,ik3,ia3
wwm_mod_xnl4v5.F90:6217:  write(luq_int,'(a1,2i3.3,a1,2i3.3,a1)') '(',ik1,ia1,'-',ik3,ia3,')'
wwm_mod_xnl4v5.F90:6218:  write(luq_int,'(4f10.4)') q_k(ik1),q_ad(ia1),q_k(ik3),q_ad(ia3)
wwm_mod_xnl4v5.F90:6219:  write(luq_int,'(4f10.4)') q_k(ik1)*cos(q_a(ia1)),q_k(ik1)*sin(q_a(ia1)),&
wwm_mod_xnl4v5.F90:6221:  write(luq_int,'(2i4)') nlocusx,19
wwm_mod_xnl4v5.F90:6353:    write(luq_int,'(1x,i4,5f8.3,11e12.4,2f11.6)')&
wwm_mod_xnl4v5.F90:6362:!!/T   write(luq_tst,'(a)') 'Q_T13V4: NSPEC'
wwm_mod_xnl4v5.F90:6364:!!/T    write(luq_tst,'(100e12.4)') (nspec(ikq,iaq),iaq=1,naq)
wwm_mod_xnl4v5.F90:6368:&  write(luq_int,'(a,4i3,i5,e13.5)') 'Q_T13V4:',ik1,ia1,ik3,ia3,nlocusx, &
wwm_mod_xnl4v5.F90:6371:!!/T if(iq_integ==3) write(luq_int,'(4i3,i5,1000e13.5)') ik1,ia1,ik3,ia3,nloc, &
wwm_mod_xnl4v5.F90:6517:!!/T    if(iq_test>=3) write(luq_tst,'(a,4f10.5)') 'Q_WEIGHT: wlog gg wlin2:', &
wwm_mod_xnl4v5.F90:6524:!!/T    if(iq_test>=3) write(luq_tst,'(a,4f10.5)') 'Q_WEIGHT: wlog gg wlin4:', &
wwm_mod_xnl4v5.F90:6576:  if(itest >= 1) write(luq_tst,'(a,i3,10f12.4)') 'Q_WEIGHT: i k2m k2a wk2 wa2 wt2(+4):',&
wwm_mod_xnl4v5.F90:6884:if(iq_screen >= 1) write(iscreen,'(a)') 'Q_XNL4V4: Checking interaction grid '
wwm_mod_xnl4v5.F90:6885:if(iq_screen >= 1) write(iscreen,'(a,2i4)') 'Q_XNL4V4: iq_search iq_type:',iq_search,iq_type
wwm_mod_xnl4v5.F90:6916:  if(iq_screen>0) write(iscreen,'(a,f12.2)') 'Q_XNL4V4: Q_SEARCHGRID called with q_depth: ',q_depth
wwm_mod_xnl4v5.F90:6919:  if(iq_screen>0) write(iscreen,'(a,f12.2)') 'Q_XNL4V4: Q_SEARCHGRID exited with q_depth: ',q_depth
wwm_mod_xnl4v5.F90:6941:  write(luq_tst,'(a)') 'NSPEC'
wwm_mod_xnl4v5.F90:6943:    write(luq_tst,'(100e12.4)') (nspec(ikq,iaq),iaq=1,naq)
wwm_mod_xnl4v5.F90:6949:    write(luq_int,'(2i4)') nkq,naq
wwm_mod_xnl4v5.F90:6950:    write(luq_int,'(10f10.3)') (q_k(ik1),ik1=1,nkq)
wwm_mod_xnl4v5.F90:6951:    write(luq_int,'(10f10.3)') (q_a(ia1)*rade,ia1=1,naq)
wwm_mod_xnl4v5.F90:6961:  if(iq_screen >= 1) write(iscreen,'(a,2i4,e12.3)') 'Q_XNL4V4: k1 nk depth:',ik1,nkq,q_depth
wwm_mod_xnl4v5.F90:6973:        if(iq_screen>=3) write(iscreen,'(a,4i4)') 'Q_XNL4V4: ik1 ia1 ik3 ia3:',ik1,ia1,ik3,ia3
wwm_mod_xnl4v5.F90:7034:!          write(*,'(a,6e13.4)') 'Check of diagonal term:',dq1,diagk1_0,t13_1,dq3,diagk3_0,t13_3
wwm_mod_xnl4v5.F90:7050:!!/R           write(iscreen,*) 'CHECK T13 QT13:',t13,qt13
wwm_mod_xnl4v5.F90:7061:            write(luq_t13,'(4i6,e13.5,2i6,f10.4)') ik1,idir1,ik3,idir3,t13,icx1,icx3,xt13
wwm_mod_xnl4v5.F90:7094:!!/F         write(luq_fil,'(a,4i3,3e11.3,2f7.2,4i2)') &
wwm_mod_xnl4v5.F90:7132:if(iq_screen>=2) write(iscreen,'(a)') 'Q_XNL4V4: Main computation ended'
wwm_netcdf.F90:154:            WRITE(WINDBG%FHNDL,*) 'YnameSec=', TRIM(Ynamesec)
wwm_netcdf.F90:155:            WRITE(WINDBG%FHNDL,*) 'lenSec=', lenSec
wwm_netcdf.F90:161:      WRITE(WINDBG%FHNDL,*) 'eStrTime=', eStrTime
wwm_netcdf.F90:398:          WRITE(DBG%FHNDL,*) 'Inconsistency in the output'
wwm_netcdf.F90:399:          WRITE(DBG%FHNDL,*) '  Vertex ', IP, ' is ', NBneighbor(IP), ' times neighbor'
wwm_netcdf.F90:403:          WRITE(DBG%FHNDL,*) 'Inconsistency in the output'
wwm_netcdf.F90:404:          WRITE(DBG%FHNDL,*) '  Vertex ', IP, ' is a neighbor'
wwm_netcdf.F90:405:          WRITE(DBG%FHNDL,*) '  but has no neighbor!'
wwm_netcdf.F90:409:          WRITE(DBG%FHNDL,*) 'Inconsistency in the output'
wwm_netcdf.F90:410:          WRITE(DBG%FHNDL,*) '  Vertex ', IP, ' has a neighbor'
wwm_netcdf.F90:411:          WRITE(DBG%FHNDL,*) '  but is not a neighbor!'
wwm_netcdf.F90:416:        WRITE(DBG%FHNDL,*) 'Find some errors in the output'
wwm_netcdf.F90:417:        WRITE(DBG%FHNDL,*) 'Please check for node contained in several boundaries'
wwm_netcdf.F90:580:        WRITE(wwmerr,*) TRIM(CallFct), ' -', idx, '-', CHRERR
wwm_nlweigt.F90:261:      WRITE(IU06,'(1H1,'' NON LINEAR INTERACTION PARAMETERS:'')')
wwm_nlweigt.F90:262:      WRITE(IU06,'(1H0,'' COMMON INDNL: CONSTANTS'')')
wwm_nlweigt.F90:263:      WRITE(IU06,*)'    ALAMD = ', ALAMD
wwm_nlweigt.F90:264:      WRITE(IU06,*)'      CON = ', CON
wwm_nlweigt.F90:265:      WRITE(IU06,*)'  DELPHI1 = ',DELPHI1
wwm_nlweigt.F90:266:      WRITE(IU06,*)'  DELPHI2 = ',DELPHI2
wwm_nlweigt.F90:267:      WRITE(IU06,'(1X,''    ACL1       ACL2   '', &
wwm_nlweigt.F90:270:      WRITE(IU06,'(1X,6F11.8)') ACL1, ACL2, CL11, CL21, DAL1, DAL2
wwm_nlweigt.F90:272:      WRITE(IU06,'(1H0,'' COMMON INDNL: FREQUENCY ARRAYS'')')
wwm_nlweigt.F90:273:      WRITE(IU06,'(1X,'' M   IKP IKP1  IKM IKM1'', &
wwm_nlweigt.F90:277:        WRITE(IU06,'(1X,I2,4I5,4F11.8,E11.3)') &
wwm_nlweigt.F90:282:      WRITE(IU06,'(1H0,'' COMMON INDNL: ANGULAR ARRAYS'')')
wwm_nlweigt.F90:283:      WRITE(IU06,'(1X,''  |--------KH = 1----------|'', &
wwm_nlweigt.F90:285:      WRITE(IU06,'(1X,'' K   K1W   K2W  K11W  K21W'', &
wwm_nlweigt.F90:288:        WRITE(IU06,'(1X,I2,8I6)') K,(K1W(K,KH), K2W(K,KH), K11W(K,KH), &
wwm_nlweigt.F90:291:      WRITE(IU06,'(1H0,'' COMMON INDNL: TAIL ARRAY FRH'')')
wwm_nlweigt.F90:292:      WRITE(IU06,'(1X,8F10.7)') (FRH(M),M=1,KFRH)
wwm_numsigma.F90:235:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_FREQUENCY'
wwm_numsigma.F90:242:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_FREQUENCY'
wwm_numtheta.F90:112:             !write(DBG%FHNDL,*) ip, is, iter, dt4di, sum(cads),sum(acq)
wwm_numtheta.F90:671:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_DIRECTION'
wwm_numtheta.F90:687:        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_DIRECTION'
wwm_output.F90:8:        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4,L5)') 'WRITING OUTPUT INTERNAL TIME', RTIME, MAIN%TMJD, OUT_HISTORY%TMJD-1.E-8, OUT_HISTORY%EMJD, (MAIN%TMJD .GE. OUT_HISTORY%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. OUT_HISTORY%EMJD)
wwm_output.F90:10:          WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'WRITING OUTPUT INTERNAL TIME', RTIME, MAIN%TMJD, OUT_HISTORY%TMJD-1.E-8, OUT_HISTORY%EMJD
wwm_output.F90:16:          WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)')  'WRITING OUTPUT INTERNAL TIME', RTIME, MAIN%TMJD, OUT_STATION%TMJD-1.E-8, OUT_STATION%EMJD
wwm_output.F90:20:        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH OUTPUT_HISTORY_AND_STATION' 
wwm_output.F90:36:        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH WWM OUTPUT'
wwm_output.F90:63:               WRITE(DBG%FHNDL,*) 'IOUTP=', VAROUT_HISTORY%IOUTP
wwm_output.F90:67:         WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH OUTPUT_HISTORY'
wwm_output.F90:82:          WRITE(STAT%FHNDL,*) 'WRITING STATION OUTPUT'
wwm_output.F90:98:        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH OUTPUT_HISTORY STATION'        
wwm_output.F90:216:           WRITE(OUT%FHNDL+1)  SNGL(TIME)
wwm_output.F90:217:           !WRITE(OUT%FHNDL+1)  (SNGL(FORCE_GLOBAL(IP,1)), SNGL(FORCE_GLOBAL(IP,2)), SNGL(OUTT_GLOBAL(IP,1))  , IP = 1, NP_GLOBAL)
wwm_output.F90:218:           WRITE(OUT%FHNDL+1)  (SNGL(OUTT_GLOBAL(IP,7)), SNGL(OUTT_GLOBAL(IP,8)), SNGL(OUTT_GLOBAL(IP,1))  , IP = 1, NP_GLOBAL)
wwm_output.F90:220:           WRITE(OUT%FHNDL+2)  SNGL(TIME)
wwm_output.F90:221:           !WRITE(OUT%FHNDL+2)  (SNGL(CURR_GLOBAL(IP,1)), SNGL(CURR_GLOBAL(IP,2)), SNGL(ZETA_SETUP(IP)), IP = 1, NP_GLOBAL)
wwm_output.F90:222:           WRITE(OUT%FHNDL+2)  (SNGL(CURR_GLOBAL(IP,1)), SNGL(CURR_GLOBAL(IP,2)), SNGL(CURR_GLOBAL(IP,5)), IP = 1, NP_GLOBAL)
wwm_output.F90:224:           WRITE(OUT%FHNDL+3)  SNGL(TIME)
wwm_output.F90:225:           WRITE(OUT%FHNDL+3)  (SNGL(WIND_GLOBAL(IP,1)), SNGL(WIND_GLOBAL(IP,2)), SNGL(OUTT_GLOBAL(IP,10))  , IP = 1, NP_GLOBAL)
wwm_output.F90:227:           WRITE(OUT%FHNDL+4)  SNGL(TIME)
wwm_output.F90:228:           WRITE(OUT%FHNDL+4)  (SNGL(OUTT_GLOBAL(IP,1)), SNGL(OUTT_GLOBAL(IP,2)), SNGL(OUTT_GLOBAL(IP,3))  , IP = 1, NP_GLOBAL)
wwm_output.F90:230:           WRITE(OUT%FHNDL+5)  SNGL(TIME)
wwm_output.F90:231:           WRITE(OUT%FHNDL+5)  (SNGL(CURR_GLOBAL(IP,1)), SNGL(CURR_GLOBAL(IP,2)), SNGL(CURR_GLOBAL(IP,5))  , IP = 1, NP_GLOBAL)
wwm_output.F90:233:           WRITE(OUT%FHNDL+6)  SNGL(TIME)
wwm_output.F90:234:           WRITE(OUT%FHNDL+6)  (SNGL(WIND_GLOBAL(IP,9)), SNGL(WIND_GLOBAL(IP,8)), SNGL(WIND_GLOBAL(IP,7))  , IP = 1, NP_GLOBAL)
wwm_output.F90:236:           WRITE(OUT%FHNDL+7)  SNGL(TIME) 
wwm_output.F90:237:           WRITE(OUT%FHNDL+7)  (SNGL(WIND_GLOBAL(IP,4)), SNGL(WIND_GLOBAL(IP,5)), SNGL(WIND_GLOBAL(IP,6))  , IP = 1, NP_GLOBAL)
wwm_output.F90:240:             WRITE(OUT%FHNDL+8) SNGL(TIME) 
wwm_output.F90:241:             WRITE(OUT%FHNDL+8)  (SNGL(ITER_GLOBAL(IP)), SNGL(ITER_GLOBAL(IP)), SNGL(ITER_GLOBAL(IP))  , IP = 1, NP_GLOBAL)
wwm_output.F90:246:               WRITE(OUT%FHNDL+9,'(10F15.6)') SNGL(WIND_GLOBAL(IP,:))
wwm_output.F90:251:             WRITE(OUT%FHNDL+10)  SNGL(TIME)
wwm_output.F90:252:             WRITE(OUT%FHNDL+10)  (SNGL(CFLCXY(IP,1)), SNGL(CFLCXY(IP,2)), SNGL(CFLCXY(IP,3)), IP = 1, MNP)
wwm_output.F90:296:         WRITE(OUT%FHNDL+1) SNGL(TIME) 
wwm_output.F90:297:         WRITE(OUT%FHNDL+1)  (SNGL(OUTT(IP,7)), SNGL(OUTT(IP,8)), SNGL(OUTT(IP,1)), IP = 1, MNP)
wwm_output.F90:299:         WRITE(OUT%FHNDL+2) SNGL(TIME)
wwm_output.F90:300:         WRITE(OUT%FHNDL+2)  (SNGL(CURR(IP,1)), SNGL(CURR(IP,2)), SNGL(DEP(IP)), IP = 1, MNP)
wwm_output.F90:302:         WRITE(OUT%FHNDL+3) SNGL(TIME)
wwm_output.F90:303:         WRITE(OUT%FHNDL+3)  (SNGL(OUTT(IP,7)), SNGL(OUTT(IP,8)), SNGL(WIND(IP,3)), IP = 1, MNP)
wwm_output.F90:305:         WRITE(OUT%FHNDL+4) SNGL(TIME)
wwm_output.F90:306:         WRITE(OUT%FHNDL+4)  (SNGL(UFRIC(IP)), SNGL(Z0(IP)), SNGL(ALPHA_CH(IP)), IP = 1, MNP)
wwm_output.F90:308:         WRITE(OUT%FHNDL+5) SNGL(TIME)
wwm_output.F90:309:         WRITE(OUT%FHNDL+5)  (SNGL(OUTT(IP,1)), SNGL(OUTT(IP,2)), SNGL(OUTT(IP,3)), IP = 1, MNP)
wwm_output.F90:311:         WRITE(OUT%FHNDL+6)  SNGL(TIME)
wwm_output.F90:312:         WRITE(OUT%FHNDL+6)  (SNGL(OUTT(IP,1)), SNGL(OUTT(IP,2)), SNGL(ISHALLOW(IP)), IP = 1, MNP)
wwm_output.F90:314:         WRITE(OUT%FHNDL+7)  SNGL(TIME)
wwm_output.F90:315:         WRITE(OUT%FHNDL+7)  (SNGL(WIND(IP,8)), SNGL(WIND(IP,9)), SNGL(WIND(IP,8)), IP = 1, MNP)
wwm_output.F90:318:           WRITE(OUT%FHNDL+8)  SNGL(TIME) 
wwm_output.F90:319:           WRITE(OUT%FHNDL+8)  (SNGL(IP_IS_STEADY(IP)), SNGL(IP_IS_STEADY(IP)), SNGL(IP_IS_STEADY(IP))  , IP = 1, NP_TOTAL)
wwm_output.F90:324:             WRITE(OUT%FHNDL+9,'(10F15.6)') SNGL(WIND(IP,:))
wwm_output.F90:329:           WRITE(OUT%FHNDL+10)  SNGL(TIME)
wwm_output.F90:330:           WRITE(OUT%FHNDL+10)  (SNGL(CFLCXY(1,IP)), SNGL(CFLCXY(2,IP)), SNGL(CFLCXY(3,IP)), IP = 1, MNP)
wwm_output.F90:335:        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH XFN_HISTORY'
wwm_output.F90:368:        WRITE(DBG%FHNDL,*) 'IP=', IP, 'HS=', HS, 'wi=', WI(I)
wwm_output.F90:374:      WRITE(DBG%FHNDL,*) 'HSinterp=', HSinterp, ' HS(b)=', HSinterpB
wwm_output.F90:375:      WRITE(DBG%FHNDL,*) 'sumAC=', sumAC
wwm_output.F90:406:      WRITE(CHRTMP,'(I2)') OUTVARS
wwm_output.F90:418:!          WRITE(DBG%FHNDL,*) 'HS=', STATION(I)%OUTPAR_NODE(1), 'sumAC=', sum(ACLOC)
wwm_output.F90:491:!         WRITE(STAT%FHNDL,*) 'DEPTH OF THE FOUND STATIONS', DEPLOC_STATIONS
wwm_output.F90:535:            WRITE(DBG%FHNDL,*) 'STATION OUT OF MESH', I
wwm_output.F90:550:            !WRITE(STAT%FHNDL,*) 'SUM BEFORE INTPAR', SUM(ACLOC_STATIONS(I,:,:))
wwm_output.F90:571:!          WRITE(STAT%FHNDL,*) I, 'SUM ACLOC 1', SUM(ACLOC)
wwm_output.F90:579:            WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
wwm_output.F90:586:          WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE(1:OUTVARS)
wwm_output.F90:599:              WRITE(OUTSP1D%FHNDL,*) MSC
wwm_output.F90:600:              WRITE(OUTSP1D%FHNDL,*) MDC
wwm_output.F90:601:              WRITE(OUTSP1D%FHNDL,*) SPSIG
wwm_output.F90:602:              WRITE(OUTSP1D%FHNDL,*) SPDIR
wwm_output.F90:604:              WRITE(OUTSP1D%FHNDL,*) STATION(I)%ISUM
wwm_output.F90:606:              WRITE(OUTSP1D%FHNDL,*) STATION(I)%IFOUND
wwm_output.F90:614:            WRITE(OUTSP1D%FHNDL,*) CTIME
wwm_output.F90:615:            WRITE(OUTSP1D%FHNDL,*) DEPLOC
wwm_output.F90:616:            WRITE(OUTSP1D%FHNDL,*) CURTXYLOC
wwm_output.F90:618:              WRITE(OUTSP1D%FHNDL,'(F15.8,3F20.10)') SPSIG(IS)/PI2,  ACOUT_1D(IS,1), ACOUT_1D(IS,2), ACOUT_1D(IS,3)
wwm_output.F90:629:              WRITE(OUTSP2D%FHNDL) MSC
wwm_output.F90:630:              WRITE(OUTSP2D%FHNDL) MDC
wwm_output.F90:631:              WRITE(OUTSP2D%FHNDL) SPSIG
wwm_output.F90:632:              WRITE(OUTSP2D%FHNDL) SPDIR
wwm_output.F90:634:              WRITE(OUTSP2D%FHNDL) STATION(I)%ISUM
wwm_output.F90:636:              WRITE(OUTSP2D%FHNDL) STATION(I)%IFOUND
wwm_output.F90:644:            WRITE(OUTSP2D%FHNDL) CTIME
wwm_output.F90:645:            WRITE(OUTSP2D%FHNDL) DEPLOC
wwm_output.F90:646:            WRITE(OUTSP2D%FHNDL) CURTXYLOC
wwm_output.F90:647:            WRITE(OUTSP2D%FHNDL) ACLOC
wwm_output.F90:648:            WRITE(OUTSP2D%FHNDL) ACOUT_2D
wwm_output.F90:656:      WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH OUTPUT_STE'
wwm_output.F90:929:        !write(DBG%FHNDL,*) 'Writing netcdf station record recs_stat=',recs_stat
wwm_output.F90:966:         WRITE(CHRTMP,'(I2)') 27
wwm_output.F90:982:             WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
wwm_output.F90:984:           WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE
wwm_output.F90:997:               WRITE(OUTSP1D%FHNDL,*) MSC, MDC
wwm_output.F90:998:               WRITE(OUTSP1D%FHNDL,*) SPSIG, SPDIR, STATION(I)%IFOUND
wwm_output.F90:1000:             WRITE(OUTSP1D%FHNDL,*) CTIME, WKLOC_STATION,DEPLOC_STATION, CURTXYLOC_STATION
wwm_output.F90:1002:               WRITE(OUTSP1D%FHNDL,'(F15.8,3F20.10)') SPSIG(IS)/PI2, ACOUT_1D(IS,1), ACOUT_1D(IS,2), ACOUT_1D(IS,3)
wwm_output.F90:1013:               WRITE(OUTSP2D%FHNDL) MSC, MDC
wwm_output.F90:1014:               WRITE(OUTSP2D%FHNDL) SPSIG, SPDIR, STATION(I)%IFOUND
wwm_output.F90:1016:             WRITE(OUTSP2D%FHNDL) CTIME, WKLOC_STATION, DEPLOC_STATION, CURTXYLOC_STATION
wwm_output.F90:1017:             WRITE(OUTSP2D%FHNDL) ACLOC, ACOUT_2D
wwm_output.F90:1028:           !WRITE(DBG%FHNDL,*) 'INTERPOLATED MYRANK =', MYRANK, I, DEPLOC(I), CURTXYLOC(I,:), SUM(WKLOC(I,:)), SUM(ACLOC_STATIONS(I,:,:))
wwm_output.F90:1031:         WRITE(DBG%FHNDL,*) 'DEPTH OF THE FOUND STATIONS LINE', DEPLOC_STATIONS
wwm_output.F90:1052:               WRITE(OUTPARM%FHNDL, TITLEFORMAT) 'TIME', OUTT_VARNAMES
wwm_output.F90:1061:               WRITE(DBG%FHNDL,*) 'STATION OUT OF MESH', I
wwm_output.F90:1071:             WRITE(OUTPARM%FHNDL,OUTPUTFORMAT) CTIME, STATION(I)%OUTPAR_NODE(1:24), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:)
wwm_output.F90:1084:                 WRITE(OUTSP1D%FHNDL,*) MSC, MDC
wwm_output.F90:1085:                 WRITE(OUTSP1D%FHNDL,*) SPSIG, SPDIR, STATION(I)%ISUM
wwm_output.F90:1087:               WRITE(OUTSP1D%FHNDL,*) CTIME, WKLOC_STATIONS(I,:), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:)
wwm_output.F90:1089:                 WRITE(OUTSP1D%FHNDL,'(F15.8,3F20.10)') SPSIG(IS)/PI2, ACOUT_1D(IS,1), ACOUT_1D(IS,2), ACOUT_1D(IS,3)
wwm_output.F90:1100:                 WRITE(OUTSP2D%FHNDL) MSC, MDC
wwm_output.F90:1101:                 WRITE(OUTSP2D%FHNDL) SPSIG, SPDIR, STATION(I)%ISUM
wwm_output.F90:1103:               WRITE(OUTSP2D%FHNDL) CTIME, WKLOC_STATIONS(I,:), DEPLOC_STATIONS(I), CURTXYLOC_STATIONS(I,:)
wwm_output.F90:1104:               WRITE(OUTSP2D%FHNDL) ACLOC, ACOUT_2D
wwm_output.F90:1167:!         WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH WINDPAR'
wwm_output.F90:1197:!         write(DBG%FHNDL,'(6F15.8)') outpar(1:6)
wwm_output.F90:1240:!         WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH INTPAR'
wwm_output.F90:1266:      WRITE(DBG%FHNDL,*) 'DEBUG AvgHS=', AvgHS, ' MaxHS=', MaxHS
wwm_output.F90:1512:         IF (ISMAX .GT. MSC) WRITE(DBG%FHNDL,*) 'ERROR IN ISMAX INTPAR_LOC'
wwm_output.F90:1588:!         WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH INTSPEC'
wwm_output.F90:1616:         WRITE(MISC%FHNDL,'(A12,6A14)') 'IP','XP','YP','DEPTH','HS','TM01','DSPR'
wwm_output.F90:1656:           WRITE(OUT%FHNDL,'(4A10)') 'ID', 'X', 'Y', 'Z'
wwm_output.F90:1658:             WRITE(4001,110) IE,',',XP(INE(1,IE)),',',YP(INE(1,IE)),',',DEP(INE(1,IE))
wwm_output.F90:1659:             WRITE(4001,110) IE,',',XP(INE(2,IE)),',',YP(INE(2,IE)),',',DEP(INE(2,IE))
wwm_output.F90:1660:             WRITE(4001,110) IE,',',XP(INE(3,IE)),',',YP(INE(3,IE)),',',DEP(INE(3,IE))
wwm_output.F90:1664:         WRITE(4002,'(4A10)') 'ID', 'X', 'Y', 'HS'
wwm_output.F90:1666:           WRITE(4002,110) IE,',',XP(INE(1,IE)),',',YP(INE(1,IE)),',',OUTT(INE(1,IE),1)
wwm_output.F90:1667:           WRITE(4002,110) IE,',',XP(INE(2,IE)),',',YP(INE(2,IE)),',',OUTT(INE(2,IE),1)
wwm_output.F90:1668:           WRITE(4002,110) IE,',',XP(INE(3,IE)),',',YP(INE(3,IE)),',',OUTT(INE(3,IE),1)
wwm_output.F90:1673:        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH OUTPUT_HISTORY_SHP'
wwm_output.F90:1786:!        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH CLSPEC'
wwm_output.F90:1837:        WRITE(STAT%FHNDL,110) TRIM(eStr), MinV, MaxV, AvgV
wwm_output.F90:1843:      WRITE(STAT%FHNDL,110) TRIM(eStr), MinV, MaxV, AvgV
wwm_output.F90:1847:!        WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH HISTORY_NC_PRINTMMA'
wwm_output.F90:2064:      WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'FINISHED WITH OUTPUT_HISTORY_NC'
wwm_parall_solver.F90:194:      WRITE(740+myrank,*) 'I5B_EXCHANGE_P4D_WWM, begin, sum(AC)=', sum(AC)
wwm_parall_solver.F90:223:      WRITE(740+myrank,*) 'Total SumErr=', SumErr
wwm_parall_solver.F90:249:      WRITE(740+myrank,*) 'SumErr=', SumErr
wwm_parall_solver.F90:253:      WRITE(740+myrank,*) 'I5B_EXCHANGE_P4D_WWM, end, sum(AC)=', sum(AC)
wwm_parall_solver.F90:564:            WRITE(700+myrank,*) 'Send IP=', IP, 'IPglob=', IP_glob
wwm_parall_solver.F90:594:            WRITE(800+myrank,*) 'Recv IP=', IPmap, 'IPglob=', IP_glob
wwm_parall_solver.F90:1039:      WRITE(740+myrank,*) 'ChromaticNr=', ChromaticNr
wwm_parall_solver.F90:1041:        WRITE(740+myrank,*) 'iVert=', iVert, 'eColor=', ListColor(iVert)
wwm_parall_solver.F90:1071:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING INIT_BLOCK_FREQDIR'
wwm_parall_solver.F90:1115:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHING INIT_BLOCK_FREQDIR'
wwm_parall_solver.F90:1141:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING INIT_BLK_L2U_ARRAY'
wwm_parall_solver.F90:1245:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHING INIT_BLK_L2U_ARRAY'
wwm_parall_solver.F90:1268:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING SYMM_INIT_COLORING'
wwm_parall_solver.F90:1272:      WRITE(740+myrank,*) 'Total residual shift=', TheRes
wwm_parall_solver.F90:1276:      WRITE(740+myrank,*) 'After COLLECT_ALL_IA_JA'
wwm_parall_solver.F90:1280:      WRITE(740+myrank,*) 'After CREATE_WWM_P2D_EXCH'
wwm_parall_solver.F90:1287:      WRITE(740+myrank,*) 'After CREATE_WWM_MAT_P2D_EXCH'
wwm_parall_solver.F90:1292:      WRITE(740+myrank,*) 'After BUILD_MULTICOLORING'
wwm_parall_solver.F90:1300:      WRITE(740+myrank,*) 'Before INIT_LOW_2_UPP_ARRAYS'
wwm_parall_solver.F90:1305:      WRITE(740+myrank,*) 'Before CALL_BLOCK_FREQDIR'
wwm_parall_solver.F90:1310:      WRITE(740+myrank,*) 'Before INIT_BLK_L2U_ARRAY'
wwm_parall_solver.F90:1315:      WRITE(740+myrank,*) 'Before COLLECT_ALL_COVLOWER'
wwm_parall_solver.F90:1320:      WRITE(740+myrank,*) 'Before INIT_COVLOWER_ARRAY'
wwm_parall_solver.F90:1325:      WRITE(740+myrank,*) 'Before DETERMINE_JSTATUS_L_U'
wwm_parall_solver.F90:1333:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHED WITH SYMM_INIT_COLORING'
wwm_parall_solver.F90:1368:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING INIT_LOW_2_UPP_ARRAYS'
wwm_parall_solver.F90:1377:        WRITE(740+myrank,*) 'I=', I, 'iRank=', iRank, 'fColor=', fColor
wwm_parall_solver.F90:1408:      WRITE(740+myrank,*) 'SIC: nbLow_recv=', nbLow_recv, ' nbUpp_send=', nbUpp_send
wwm_parall_solver.F90:1464:      WRITE(740+myrank,*) 'SumErr(meth1/meth2) CovLower=', SumErr
wwm_parall_solver.F90:1465:      WRITE(740+myrank,*) 'MNP=', MNP, ' sum(CovLower)=', sum(CovLower)
wwm_parall_solver.F90:1468:          WRITE(740+myrank,*) IP, CovLower(IP), CovLower_meth2(IP)
wwm_parall_solver.F90:1472:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHED WITH INIT_LOW_2_UPP_ARRAYS'
wwm_parall_solver.F90:1503:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING INIT_COVLOWER_ARRAY'
wwm_parall_solver.F90:1527:      WRITE(740+myrank,*) 'nbMap0=', nbMap0, ' nbMap1=', nbMap1
wwm_parall_solver.F90:1545:      WRITE(740+myrank,*) 'U2L eColor=', eColor
wwm_parall_solver.F90:1551:        WRITE(740+myrank,*) 'U2L iNeigh=', iNeigh, ' iProc=', iProc, 'fColor=', fColor
wwm_parall_solver.F90:1585:          WRITE(740+myrank,*) '   U2L nbCommon_recv=', nbCommon_recv
wwm_parall_solver.F90:1609:          WRITE(740+myrank,*) '   U2L nbCommon_send=', nbCommon_send
wwm_parall_solver.F90:1614:      WRITE(740+myrank,*) 'WWM_P2D: u2l_nnbr_send=', u2l_nnbr_send
wwm_parall_solver.F90:1615:      WRITE(740+myrank,*) 'WWM_P2D: u2l_nnbr_recv=', u2l_nnbr_recv
wwm_parall_solver.F90:1645:      WRITE(740+myrank,*) 'WWM_P2D: wwm_nnbr=', wwm_nnbr
wwm_parall_solver.F90:1646:      WRITE(740+myrank,*) 'WWM_P2D: wwm_ListNeigh built'
wwm_parall_solver.F90:1651:      WRITE(740+myrank,*) 'WWM_P2D: alloc done'
wwm_parall_solver.F90:1767:      WRITE(740+myrank,*) 'wwm_nnbr=', wwm_nnbr
wwm_parall_solver.F90:1768:      WRITE(740+myrank,*) 'sum(ListCovLower)=', sum(LocalColor % ListCovLower)
wwm_parall_solver.F90:1822:        WRITE(740+myrank,*) '   nbCommon(send/recv)=', nbCommon_send, nbCommon_recv
wwm_parall_solver.F90:1826:      WRITE(740+myrank,*) 'WWM_P2D: sync_nnbr_send=', sync_nnbr_send
wwm_parall_solver.F90:1827:      WRITE(740+myrank,*) 'WWM_P2D: sync_nnbr_recv=', sync_nnbr_recv
wwm_parall_solver.F90:1855:      WRITE(740+myrank,*) 'SYNC sync_nnbr_send=', sync_nnbr_send
wwm_parall_solver.F90:1875:        WRITE(740+myrank,*) '   SYNC iNeigh=', iNeigh, ' nbCommon=', nbCommon
wwm_parall_solver.F90:1897:              WRITE(740+myrank,*) 'idx=', idx, 'IP=', IP
wwm_parall_solver.F90:1898:              WRITE(740+myrank,*) '  IP_glob=', IP_glob, 'IPmap=', IPmap
wwm_parall_solver.F90:1967:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHED WITH INIT_COVLOWER_ARRAY'
wwm_parall_solver.F90:2032:            WRITE(myrank+919,*) 'MNP=', MNP, 'NP_RES=', NP_RES
wwm_parall_solver.F90:2033:            WRITE(myrank+919,*) 'IP=', IP, ' JP=', JP
wwm_parall_solver.F90:2034:            WRITE(myrank+919,*) 'IPcovLower=', LocalColor%CovLower(IP)
wwm_parall_solver.F90:2035:            WRITE(myrank+919,*) 'JPcovLower=', LocalColor%CovLower(JP)
wwm_parall_solver.F90:2036:            WRITE(myrank+919,*) 'We have major error'
wwm_parall_solver.F90:2041:      WRITE(740+myrank,*) 'sum(Jstatus_L)=', sum(Jstatus_L)
wwm_parall_solver.F90:2042:      WRITE(740+myrank,*) 'sum(Jstatus_U)=', sum(Jstatus_U)
wwm_parall_solver.F90:2325:      !WRITE(myrank+640,*) 'MIN ASPAR_block=', minval(SolDat % ASPAR_block)
wwm_parall_solver.F90:2539:!                  WRITE(740+myrank,*) 'hVal=', hVal, 'ACtest=', ACtest(IS,ID,JP)
wwm_parall_solver.F90:2540:!                  WRITE(740+myrank,*) '      diff=', hVal - ACtest(IS,ID,JP)
wwm_parall_solver.F90:2546:!                  WRITE(740+myrank,*) '1J=', J, 'eCoeff=', eCoeff, 'eCoeffB=', eCoeffB
wwm_parall_solver.F90:2547:!                  WRITE(740+myrank,*) '      diff=', eCoeff - eCoeffB
wwm_parall_solver.F90:2621:!                  WRITE(740+myrank,*) '2J=', J, 'eCoeff=', eCoeff, 'eCoeffB=', eCoeffB
wwm_parall_solver.F90:2622:!                  WRITE(740+myrank,*) '      diff=', eCoeff - eCoeffB
wwm_parall_solver.F90:2636:!              WRITE(740+myrank,*) '3J=', J, 'eCoeff=', eCoeff, 'eCoeffB=', eCoeffB
wwm_parall_solver.F90:2637:!              WRITE(740+myrank,*) '      diff=', eCoeff - eCoeffB
wwm_parall_solver.F90:2744:      WRITE(740+myrank,*) 'ACret   NP_RES cohenrency error=', Lerror
wwm_parall_solver.F90:2746:      WRITE(740+myrank,*) 'ACtest1 NP_RES cohenrency error=', Lerror
wwm_parall_solver.F90:2748:      WRITE(740+myrank,*) 'ACtest2 NP_RES cohenrency error=', Lerror
wwm_parall_solver.F90:2751:      WRITE(740+myrank,*) 'ACret   MNP cohenrency error=', Lerror
wwm_parall_solver.F90:2753:      WRITE(740+myrank,*) 'ACtest1 MNP cohenrency error=', Lerror
wwm_parall_solver.F90:2755:      WRITE(740+myrank,*) 'ACtest2 MNP cohenrency error=', Lerror
wwm_parall_solver.F90:3225:!      WRITE(740+myrank,*) 'Beginning solution'
wwm_parall_solver.F90:3229:!        WRITE(740+myrank,*) 'nbIter=', nbIter
wwm_parall_solver.F90:3298:!        WRITE(740+myrank,*) 'CritVal=', CritVal
wwm_parall_solver.F90:3313:!      WRITE(740+myrank,*) 'End BCGS_REORG'
wwm_parall_solver.F90:3398:      WRITE(myrank+740,*) 'After SYMM_INIT_COLORING'
wwm_parall_solver.F90:3403:      WRITE(myrank+740,*) 'After I5B_ALLOCATE'
wwm_parall_solver.F90:3466:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING WWM_SOLVER_EIMPS'
wwm_parall_solver.F90:3475:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHING WWM_SOLVER_EIMPS'
wwm_parall_solver.F90:3728:      WRITE(3000+myrank,*) 'iMSCblock=', iMSCblock
wwm_parall_solver.F90:3729:      WRITE(3000+myrank,*) 'IS12=', IS1, IS2
wwm_parall_solver.F90:3731:        WRITE(3000+myrank,*) 'A: IS, sum(ASPAR)=', IS, sum(ASPAR(IS,:,:))
wwm_parall_solver.F90:3750:      WRITE(3000+myrank,*) 'iMSCblock=', iMSCblock
wwm_parall_solver.F90:3751:      WRITE(3000+myrank,*) 'IS12=', IS1, IS2
wwm_parall_solver.F90:3753:        WRITE(3000+myrank,*) 'B: IS, sum(ASPAR)=', IS, sum(ASPAR(IS,:,:))
wwm_parall_solver.F90:3768:      WRITE(3000+myrank,*) 'sum(ASPAR )=', sum(ASPAR)
wwm_parall_solver.F90:3769:      WRITE(3000+myrank,*) 'sum(B     )=', sum(B)
wwm_parall_solver.F90:3770:      WRITE(3000+myrank,*) 'iMSCblock=', iMSCblock
wwm_parall_solver.F90:3771:      WRITE(3000+myrank,*) 'IS12=', IS1, IS2
wwm_parall_solver.F90:3773:        WRITE(3000+myrank,*) 'C: IS, sum(ASPAR)=', IS, sum(ASPAR(IS,:,:))
wwm_parall_solver.F90:3792:      WRITE(740+myrank,*) 'Begin of EIMPS_B_BLOCK'
wwm_parall_solver.F90:3814:      WRITE(3000+myrank,*)  'sum(B     )=', sum(B)
wwm_parall_solver.F90:3850:      WRITE(740+myrank,*) 'Begin of EIMPS_ASPAR_B_BLOCK'
wwm_parall_solver.F90:3866:      WRITE(740+myrank,*) ' Before MNE loop'
wwm_parall_solver.F90:4048:!      WRITE(*,*) SUM(IMATRAA), SUM(IMATDAA)
wwm_parall_solver.F90:4051:      WRITE(3000+myrank,*)  'sum(ASPAR )=', sum(ASPAR)
wwm_parall_solver.F90:4052:      WRITE(3000+myrank,*)  'sum(B     )=', sum(B)
wwm_parall_solver.F90:4054:        WRITE(3000+myrank,*) 'IS, sum(ASPAR)=', IS, sum(ASPAR(IS,:,:))
wwm_parall_solver.F90:4092:      WRITE(740+myrank,*) 'Begin of EIMPS_ASPAR_B_BLOCK'
wwm_parall_solver.F90:4108:      WRITE(740+myrank,*) ' Before MNE loop'
wwm_parall_solver.F90:4248:      WRITE(3000+myrank,*)  'sum(ASPAR )=', sum(ASPAR)
wwm_parall_solver.F90:4249:      WRITE(3000+myrank,*)  'sum(B     )=', sum(B)
wwm_parall_solver.F90:4251:        WRITE(3000+myrank,*) 'IS, sum(ASPAR)=', IS, sum(ASPAR(IS,:,:))
wwm_parall_solver.F90:4274:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING I5B_EIMPS'
wwm_parall_solver.F90:4278:      WRITE(740+myrank,*) 'Begin I5B_EIMPS'
wwm_parall_solver.F90:4296:      WRITE(740+myrank,*) 'After ASPAR init'
wwm_parall_solver.F90:4306:      WRITE(740+myrank,*) 'After EXCHANGE_P4D_WWM'
wwm_parall_solver.F90:4310:      WRITE(740+myrank,*) 'After I5B_EXCHANGE_ASPAR'
wwm_parall_solver.F90:4317:      WRITE(740+myrank,*) 'After I5B_CREATE_PRECOND'
wwm_parall_solver.F90:4331:      WRITE(STAT%FHNDL,*) 'nbIter=', nbIter, 'L2/LINF=', maxval(Norm_L2), maxval(Norm_LINF)
wwm_parall_solver.F90:4334:        Write(myrank+591,*) 'Clearing ENDING'
wwm_parall_solver.F90:4338:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHING I5B_EIMPS'
wwm_parall_solver.F90:4364:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'ENTERING I4_EIMPS'
wwm_parall_solver.F90:4369:        WRITE(240+myrank,*) 'iMSCblock=', iMSCblock
wwm_parall_solver.F90:4404:      WRITE(STAT%FHNDL,*) 'nbIter=', nbIterMax, 'L2/LINF=', Max_L2, Max_LINF
wwm_parall_solver.F90:4405:      WRITE(STAT%FHNDL,'("+TRACE......",A)') 'FINISHED I4_EIMPS'
wwm_parall_solver.F90:4454:!      WRITE(*,*) SUM(IMATRAA), SUM(IMATDAA)
wwm_parall_solver.F90:4611:            !write(*,'(I10,4F25.20,L10)') ip, p_is_converged, solverthr, sumu, sumx, p_is_converged .lt. solverthr
wwm_parall_solver.F90:4625:             WRITE(850+myrank,'(3I10,2F20.17,L10)') NBITER, IP, IPLG(IP), p_is_converged, solverthr, p_is_converged .lt. solverthr
wwm_parall_solver.F90:4636:          !if (myrank == 0) write(*,*) nbiter, is_converged, np_global, p_is_converged, solverthr
wwm_parall_solver.F90:4694:        WRITE(STAT%FHNDL,'(A10,3I10,2F20.10)') 'solver', nbiter, is_converged, np_global-is_converged, p_is_converged, pmin
wwm_parall_solver.F90:4695:        !WRITE(*,'(A10,4I10,2F20.10)') 'solver', nbiter, maxiter, is_converged, np_global-is_converged, p_is_converged, pmin
wwm_parall_solver.F90:4723:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPROCESSING SOURCES AND ADVECTION  ', TIME2-TIME1
wwm_parall_solver.F90:4724:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPROCESSING REFRACTION             ', TIME3-TIME2
wwm_parall_solver.F90:4725:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'ITERATION                            ', TIME4-TIME3
wwm_parall_solver.F90:4726:        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'STORE RESULT                         ', TIME5-TIME4
wwm_petsc_block.F90:236:!        write(DBG%FHNDL,*) "openmp threads", nOMPthreads
wwm_petsc_block.F90:462:              WRITE(740+myrank,*) 'IS/ID/IP=', ISS, IDD, IPpetsc
wwm_petsc_block.F90:463:              WRITE(740+myrank,*) 'idxpos=', idxpos
wwm_petsc_block.F90:626:!         write(DBG%FHNDL,*) rank, "petsc_small nnz_new", nnz_new, " old", NNZ
wwm_petsc_block.F90:627:!         write(DBG%FHNDL,*) rank, "o_nnz_new", o_nnz_new
wwm_petsc_block.F90:1090:            WRITE(640+myrank,*) 'idxposfinal=', idxpos
wwm_petsc_block.F90:1288:                !if (value(1)/max(thr,value(2)) .lt. 100. .or. value(1)/max(thr,value(3)) .lt. 100.) write(*,*) value(1), value(1)/max(thr,value(2)), value(1)/max(thr,value(3))
wwm_petsc_block.F90:1310:!            WRITE(640+myrank,*) 'IS/ID/IP=', ISS,IDD, IPpetsc, idxpos
wwm_petsc_block.F90:1527:          !WRITE(DBG%FHNDL,*) IP, DIFRM(IP), C(:,IP)
wwm_petsc_block.F90:1660:          !WRITE(DBG%FHNDL,*) IP, DIFRM(IP), C(:,IP)
wwm_petsc_block.F90:1739:          !write(stat%fhndl,*) 'Failure to converge'
wwm_petsc_block.F90:1745:             write(DBG%FHNDL,*) "Failure to converge\n"
wwm_petsc_block.F90:1749:            if(iteration /= 0)  write(DBG%FHNDL,*) "Number of iterations", iteration
wwm_petsc_block.F90:1895:          write(DBG%FHNDL,*) rank, "aspar2petscAspar < 1 !! IPpetsclocal IS ID asparposi", IPpetscLocal, ISS, IDD, asparPosition
wwm_petsc_block.F90:1934:          write(DBG%FHNDL,*) rank, "oAspar2petscAspar < 1 !! IPpetsclocal IS ID asparposi", IPpetscLocal, ISS, IDD, asparPosition
wwm_petsc_block.F90:2093:          write(DBG%FHNDL,*) "check matrix diagonal Accuracy"
wwm_petsc_block.F90:2094:          if(present(ISS) .and. present(IDD)) write(DBG%FHNDL,*) "ISS IDD", ISS, IDD
wwm_petsc_block.F90:2095:          write(DBG%FHNDL,*) "minimum at (big matrix row)" , positionMin, ": ", valueMin
wwm_petsc_block.F90:2096:          write(DBG%FHNDL,*) "maximum at (big matrix row)" , positionMax, ": ", valueMax
wwm_petsc_block.F90:2097:          write(DBG%FHNDL,*) "mean" , summe / globalSize
wwm_petsc_block.F90:2099:          write(DBG%FHNDL,*) "first 10 entries which are smaller than", epsilon
wwm_petsc_block.F90:2100:          write(DBG%FHNDL,*) "bigmatrix | IPpetsc global | APP global |  IOBP | IOBPD    |    ISS    |    IDD"
wwm_petsc_block.F90:2102:            write(DBG%FHNDL,*) entriesDetail(i,:)
wwm_petsc_block.F90:2105:          write(DBG%FHNDL,*) rank, " There are total ", zeroElementsCounter," entries"
wwm_petsc_block.F90:2107:          write(DBG%FHNDL,*) "check matrix diagonal Accuracy Ende. Time: ", endTime - startTime," sec"
wwm_petsc_controller.F90:37:          write(DBG%FHNDL,*) "PETSC_INIT() you can use AMETHOD=4 "
wwm_petsc_controller.F90:38:          write(DBG%FHNDL,*) "for petsc parallel or AMETHOD=5 for petsc block."
wwm_petsc_controller.F90:39:          write(DBG%FHNDL,*) "Other AMETHOD numbers are for other solvers."
wwm_petsc_controller.F90:40:          write(DBG%FHNDL,*) "AMETHOD =" ,  AMETHOD
wwm_petsc_controller.F90:92:            write(DBG%FHNDL,*) "PETSC_FINALIZE() you can use AMETHOD=4"
wwm_petsc_controller.F90:93:            write(DBG%FHNDL,*) "for petsc parallel or AMETHOD=5 for petsc block."
wwm_petsc_controller.F90:94:            write(DBG%FHNDL,*) "Other AMETHOD numbers are for other solvers."
wwm_petsc_controller.F90:95:            write(DBG%FHNDL,*) "AMETHOD =" , AMETHOD
wwm_petsc_parallel.F90:164:!         write(DBG%FHNDL,*) rank, "nnz_new", nnz_new, " old", NNZ
wwm_petsc_parallel.F90:165:!         write(DBG%FHNDL,*) rank, "o_nnz_new", o_nnz_new
wwm_petsc_parallel.F90:180:          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
wwm_petsc_parallel.F90:186:          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
wwm_petsc_parallel.F90:257:          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
wwm_petsc_parallel.F90:521:           !write(stat%fhndl,*) 'Failure to converge'
wwm_petsc_parallel.F90:527:              write(DBG%FHNDL,*) "Failure to converge\n"
wwm_petsc_parallel.F90:535:               write(DBG%FHNDL,*) "mean number of iterations", iterationSum / real((MSC*MDC))
wwm_petscpool.F90:239:          write(DBG%FHNDL,*) "check ASPAR diagonal Accuracy"
wwm_petscpool.F90:240:          if(present(ISS) .and. present(IDD)) write(DBG%FHNDL,*) "ISS IDD", ISS, IDD
wwm_petscpool.F90:241:          write(DBG%FHNDL,*) "minimum at (IP global)" , iplg(positionMin), ": ", valueMin
wwm_petscpool.F90:242:          write(DBG%FHNDL,*) "maximum at (IP global)" , iplg(positionMax), ": ", valueMax
wwm_petscpool.F90:243:          write(DBG%FHNDL,*) "mean" , mean
wwm_petscpool.F90:245:          write(DBG%FHNDL,*) "first 10 entries which are smaller than", epsilon
wwm_petscpool.F90:246:          write(DBG%FHNDL,*) "IP (global)\tIOBP"
wwm_petscpool.F90:248:            write(DBG%FHNDL,*) entriesDetail(i,:)
wwm_petscpool.F90:251:          write(DBG%FHNDL,*) rank, " There are total ",zeroElementsCounter," entries"
wwm_petscpool.F90:253:          write(DBG%FHNDL,*) "check ASPAR diagonal Accuracy Ende. Time: ",endTime - startTime," sec"
wwm_petscpool.F90:333:        write(DBG%FHNDL,*) "relative convergence tolerance", rtol
wwm_petscpool.F90:334:        write(DBG%FHNDL,*) "absolute convergence tolerance", abstol
wwm_petscpool.F90:335:        write(DBG%FHNDL,*) "divergence tolerance", dtol
wwm_petscpool.F90:336:        write(DBG%FHNDL,*) "maximum number of iterations", maxits
wwm_petscpool.F90:358:        write(DBG%FHNDL,*) "using KSP: ", trim(ksp)
wwm_petscpool.F90:362:          write(DBG%FHNDL,*) "KSP using SAME_PRECONDITIONER"
wwm_petscpool.F90:364:          write(DBG%FHNDL,*) "KSP using SAME_NONZERO_PATTERN"
wwm_petscpool.F90:366:          write(DBG%FHNDL,*) "KSP using DIFFERENT_NONZERO_PATTERN"
wwm_petscpool.F90:368:          write(DBG%FHNDL,*) "KSP using SUBSET_NONZERO_PATTERN"
wwm_petscpool.F90:370:          write(DBG%FHNDL,*) "KSPGetOperators() unknown operator set!"
wwm_petscpool.F90:375:        write(DBG%FHNDL,*) "using PC: ", trim(pc)
wwm_petscpool.F90:406:            if(rank == 0) write(DBG%FHNDL,*) "IS ID", ISS, IDD
wwm_petscpool.F90:411:            if(rank == 0) write(DBG%FHNDL,*) "matrix:"
wwm_petscpool.F90:417:            if(rank == 0) write(DBG%FHNDL,*) "X"
wwm_petscpool.F90:423:            if(rank == 0) write(DBG%FHNDL,*) "rhs"
wwm_petscpool.F90:439:        write(DBG%FHNDL,*) "block size", matInfo(MAT_INFO_BLOCK_SIZE)
wwm_petscpool.F90:440:        write(DBG%FHNDL,*) "number of nonzeros allocated", matInfo(MAT_INFO_NZ_ALLOCATED)
wwm_petscpool.F90:441:        write(DBG%FHNDL,*) "number of nonzeros used", matInfo(MAT_INFO_NZ_USED)
wwm_petscpool.F90:442:        write(DBG%FHNDL,*) "number of nonzeros uneeded", matInfo(MAT_INFO_NZ_UNNEEDED)
wwm_petscpool.F90:443:        write(DBG%FHNDL,*) "memory allocated", matInfo(MAT_INFO_MEMORY)
wwm_petscpool.F90:444:        write(DBG%FHNDL,*) "number of matrix assemblies called ", matInfo(MAT_INFO_ASSEMBLIES)
wwm_petscpool.F90:445:        write(DBG%FHNDL,*) "number of mallocs during MatSetValues()", matInfo(MAT_INFO_MALLOCS)
wwm_petscpool.F90:446:        write(DBG%FHNDL,*) "fill ratio for LU/ILU given", matInfo(MAT_INFO_FILL_RATIO_GIVEN)
wwm_petscpool.F90:447:        write(DBG%FHNDL,*) "fill ratio for LU/ILU needed", matInfo(MAT_INFO_FILL_RATIO_NEEDED)
wwm_petscpool.F90:448:        write(DBG%FHNDL,*) "number of mallocs during factorization", matInfo(MAT_INFO_FACTOR_MALLOCS)
wwm_petscpool.F90:476:        write(DBG%FHNDL,*) "global matrix properties"
wwm_petscpool.F90:477:        write(DBG%FHNDL,*) "NORM 1/2/inf"
wwm_petscpool.F90:478:        write(DBG%FHNDL,*) norm1, norm2, norminf
wwm_petscpool.F90:479:        write(DBG%FHNDL,*) "diagMin, diagMax, ratio"
wwm_petscpool.F90:480:        write(DBG%FHNDL,*) diagMin, diagMax, diagMax/diagMin
wwm_petscpool.F90:481:        write(DBG%FHNDL,*) "local matrix properties"
wwm_petscpool.F90:500:      write(DBG%FHNDL,*) rank, "NORM1", norm1
wwm_petscpool.F90:501:      write(DBG%FHNDL,*) rank, "NORM2", norm2
wwm_petscpool.F90:502:      write(DBG%FHNDL,*) rank, "NORMinf", norminf
wwm_petscpool.F90:503:      write(DBG%FHNDL,*) rank, "diagMin", diagMin
wwm_petscpool.F90:504:      write(DBG%FHNDL,*) rank, "diagMax", diagMax
wwm_petscpool.F90:505:      write(DBG%FHNDL,*) rank, "diagRatio", diagMax/diagMin
wwm_petscpool.F90:578:          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
wwm_petscpool.F90:605:          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
wwm_petscpool.F90:621:          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
wwm_petscpool.F90:641:          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
wwm_petscpool.F90:660:          write(DBG%FHNDL,*) __FILE__, " Line", __LINE__
wwm_petscpool.F90:681:!           write(DBG%FHNDL,*) rank, "Global Number of Nodes" , np_global
wwm_petscpool.F90:682:!           write(DBG%FHNDL,*) rank, "Local Number of resident nodes", NP_RES
wwm_petscpool.F90:683:!           write(DBG%FHNDL,*) rank, "Local Number of ghost nodes", npg
wwm_petscpool.F90:684:!           write(DBG%FHNDL,*) rank, "local Number of nodes in augmented subdomain (NP_RES+npg)", MNP
wwm_petscpool.F90:685:!           write(DBG%FHNDL,*) rank, "Local Number of nodes without interface and ghost nodes", nNodesWithoutInterfaceGhosts
wwm_petscpool.F90:686:!           write(DBG%FHNDL,*) rank, "Local Number of ghost + interface nodes", nghost
wwm_petscpool.F90:741:          WRITE(CHK%FHNDL, NML=PETScOptions)
wwm_petscpool.F90:755:          write(DBG%FHNDL,*) "Strange input in namelist PETScOptions"
wwm_petscpool.F90:758:        if(rtolStrage)    write(DBG%FHNDL,*) "rtol > 1. Are you sure?", rtol
wwm_petscpool.F90:759:        if(abstolStrange) write(DBG%FHNDL,*) "abstol > 1. Are you sure?", abstol
wwm_petscpool.F90:760:        if(dtolStrange)   write(DBG%FHNDL,*) "dtol < 1. Are you sure?", dtol
wwm_petscpool.F90:761:        if(maxitsStrange) write(DBG%FHNDL,*) "maxits < 1 ir maxits > 10000 Are you sure?", maxits
wwm_petscpool.F90:765:!         if(rank == 0) write(CHK%FHNDL, NML=PETScOptions)
wwm_petscpool.F90:781:           write(DBG%FHNDL,*) __FILE__, " Line ", __LINE__ ," petsc matrix was not created. call createMatrix() befor createSolver()"
wwm_petsc_seriell.F90:264:         !WRITE(*,*) TIME2-TIME1, TIME3-TIME2, TIME4-TIME3, TIME5-TIME4
wwm_petsc_seriell.F90:268:           write(*,*) "Failure to converge\n"
wwm_petsc_seriell.F90:272:            if(iterationen /= 0)  write(*,*) "Number of iterations", iss,idd,iterationen
wwm_petsc_seriell.F90:273:!            write(*,*) "Number of iterations", iterationen
wwm_sdiss_ardh_vec.F90:485:            WRITE(111116,'(10F15.8)') D(IJ,K,M),DTEMP,SIG(M),RENEWALFREQ(IJ,K),BTH0(IJ,M),BTH(IJ,K,M)
wwm_sdissip.F90:104:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111119,*) '------- STARTING DISSIPATION -------'
wwm_sdissip.F90:114:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111119,'(5F20.10)') SUM(F), SUM(FL), SUM(SL), &
wwm_sdissip.F90:159:          WRITE(IU06,*) '   SUB. SDISSIP: START DO-LOOP (ISHALLO=0)'
wwm_sdissip.F90:190:          IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111119,'(5F20.10)')SDS(IJ),TEMP1(IJ),&
wwm_sdissip.F90:212:               !IF (IJ == 339) write(*,'(I10,8F20.10)') IJ, DEP(IJ), SDS(IJ), Q, COEF_B_J, ALPH, F1MEAN(IJ), 4*SQRT(EMAX), 4*SQRT(EMEAN(IJ))
wwm_sdissip.F90:233:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111119,'(2F30.20)') SUM(FL), SUM(SL)
wwm_sdissip.F90:234:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111119,*) '------- FINISHED DISSIPATION -------' 
wwm_serv_xnl4v5.F90:223:      write(stat%fhndl,*) 'z_intp1: i1 x1(i1) x1(i1+1):',i1,x1(i1),x1(i1+1)
wwm_serv_xnl4v5.F90:231:      write(stat%fhndl,*) 'z_intp1: i2 x2(i2) x2(i2+1):',i2,x2(i2),x2(i2+1)
wwm_serv_xnl4v5.F90:240:      write(stat%fhndl,*) 'z_intp1: i1 x1(i1) x1(i1+1):',i1,x1(i1),x1(i1+1)
wwm_serv_xnl4v5.F90:248:      write(stat%fhndl,*) 'z_intp1: i2 x2(i2) x2(i2+1):',i2,x2(i2),x2(i2+1)
wwm_serv_xnl4v5.F90:540:    write(*,'(a,i4)') 'Z_ROOT2: invalid unit number:',iprint
wwm_serv_xnl4v5.F90:555:!if(luprint > 0) write(luprint,'(a,4e13.5)') &
wwm_serv_xnl4v5.F90:570:!      if(luprint>0) write(luprint,'(a,4e13.5)') &
wwm_serv_xnl4v5.F90:574:!        if(luprint>0) write(luprint,'(a)') 'Z_ROOT2: xnew=zriddr'
wwm_serv_xnl4v5.F90:600:      if(luprint > 0) write(luprint,'(a,i4,5e14.6)') &
wwm_serv_xnl4v5.F90:605:   if(luprint > 0) write(luprint,'(a)') 'Z_ROOT2: -> ierr=2'
wwm_serv_xnl4v5.F90:621:if(luprint > 0) write(luprint,'(a,2i3,5e13.5)') &
wwm_sinput.F90:151:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,*) '------- STARTING SINPUT --------'
wwm_sinput.F90:153:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,'(10F15.7)') CONST3,XKAPPAD,CONST3,SUM(F)
wwm_sinput.F90:168:          !IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,'(2I10,5F20.10)') M,K,TH(K),THWNEW(IJ)
wwm_sinput.F90:232:           !IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,'(5F30.20)') UCN1(IJ),UCN2(IJ),ZCN(IJ)
wwm_sinput.F90:233:           !IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,'(5F30.20)') CNSN(IJ),XV1D(IJ),XV2D(IJ)
wwm_sinput.F90:269:            !IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,'(2I10,10F15.7)') M, K, TEMPD(IJ,K) 
wwm_sinput.F90:285:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,'(A30)') '--------- NOW THE SPECTRA ---------'
wwm_sinput.F90:289:          IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,'(3F30.20)') SUM(F(IJ,:,M)), SUM(FL(IJ,:,M)), SUM(SL(IJ,:,M))
wwm_sinput.F90:293:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111114,*) '------- FINISHED SINPUT ---------'
wwm_snl3.F90:41:      !if (ip == 1786) write(stat%fhndl,'(A20,I10,8F15.10)') 'URSELL',IP,DEP(IP),HS,SMESPC,URSELL,(G9 * HS),(TWO*SQRT(TWO)*SMESPC**2*DEP(IP)**2)
wwm_snl3.F90:43:!      write(*,*) '---- calling snl3 -----', ip, iobp(ip)
wwm_snl3.F90:54:         WRITE(DBG%FHNDL,*) 'CHECK RESOLUTION', IRES, FACSCL, FACRES, SIGLOW
wwm_snl3.F90:113:              IF (IP == 1786 .AND. SA(IS,ID) .GT. THR) WRITE(*,'(2I10,10F25.10)') IS, ID, DEP(IP), URSELL, PTRIAD(5), RINT**2, SINBPH, SA(IS,ID), FT, ( EM * (EM - 2*E0 ))
wwm_snl3.F90:122:              !IF (IP == 1786)  WRITE(*,'(2I10,4F15.10,I10)') IS, ID, STRI, SA(IS,ID), SA(IS+ISP1,ID) , SA(IS+ISP,ID), ISP+IS
wwm_snl3.F90:140:        WRITE(*,*) 'FINAL SUMS', SUM(IMATRA), SUM(IMATDA), SUM(SSNL3)
wwm_snl3.F90:278:!        write(*,*) 'delta', res 
wwm_snl3.F90:299:      if (ip == 157) write(*,*) 'ddelta_dx', res
wwm_snl3.F90:358:        if (ip == 157) write(*,'(A20,7F15.10)') 'w    ', res, tau(ip, is, is1, is2,n1), a, b, c, d, e
wwm_snl3.F90:441:        !if (ip == 157) write(*,*) 'dwdx', res 
wwm_snl3.F90:457:      if (ip == 157) write(stat%fhndl,*) 'k', res
wwm_snl3.F90:474:       if (ip == 157) write(stat%fhndl,*) 'j', res 
wwm_snl3.F90:540:         if (ip == 157) write(*,'(2I10,A10,3F20.15)') is, id, 'super',fr(is),super,superd
wwm_snl3.F90:560:         if (ip == 157) write(*,'(2I10,A10,3F20.15)') is,id,'sub',fr(is),sub,subd
wwm_snl3.F90:570:         write(stat%fhndl,'(2i10,3f15.10)') ip, is, fr(is), sum(snl3(is,:)) * DDIR
wwm_snl3.F90:587:        write(stat%fhndl,*) 'kron-delta', i, j, res
wwm_snl3.F90:975:!        WRITE(*,*) DEP_2, DEP_3, I2, I1, XIS, XISLN, ISP, ISP1, WISP, WISP1, ISM, ISM1, WISM, WISM1
wwm_snl3.F90:976:!        WRITE(*,*) SUM(IMATRA), SUM(IMATDA), SUM(SSNL3)
wwm_snl3.F90:1215:      !if (ip == 1786) write(stat%fhndl,'(A20,I10,8F15.10)') 'URSELL',IP,DEP(IP),HS,SMESPC,URSELL,(G9 * HS),(TWO*SQRT(TWO)*SMESPC**2*DEP(IP)**2)
wwm_snl3.F90:1217:!      write(*,*) '---- calling snl3 -----', ip, iobp(ip)
wwm_snl3.F90:1228:         WRITE(DBG%FHNDL,*) 'CHECK RESOLUTION', IRES, FACSCL, FACRES, SIGLOW
wwm_snl3.F90:1291:              !IF (IP == 1786 .AND. SA(IS,ID) .GT. THR) WRITE(*,'(2I10,10F25.10)') IS, ID, DEP(IP), URSELL, PTRIAD(5), RINT**2, SINBPH, SA(IS,ID), FT, ( EM * (EM - 2*E0 ))
wwm_snl42.F90:141:            WRITE(STAT % FHNDL,*) 'PARAMETERS FOR SNL'
wwm_snl42.F90:142:            WRITE(STAT % FHNDL,*) IDP, IDP1, IDM, IDM1
wwm_snl42.F90:143:            WRITE(STAT % FHNDL,*) ISP, ISP1, ISM, ISM1
wwm_snl42.F90:144:            WRITE(STAT % FHNDL,*) ISLOW, ISHGH, ISCLW, ISCHG, IDLOW, IDHGH
wwm_snl42.F90:145:            WRITE(STAT % FHNDL,*) XIS
wwm_snl42.F90:146:            WRITE(STAT % FHNDL,*) AWG1, AWG2, AWG3, AWG4
wwm_snl42.F90:147:            WRITE(STAT % FHNDL,*) AWG5, AWG6, AWG7, AWG8
wwm_snl42.F90:148:            WRITE(STAT % FHNDL,*) '---------------------------------------'
wwm_snl4.F90:140:            WRITE(STAT%FHNDL,*) 'PARAMETER 4 SNL'
wwm_snl4.F90:141:            WRITE(STAT%FHNDL,*) IDP, IDP1, IDM, IDM1
wwm_snl4.F90:142:            WRITE(STAT%FHNDL,*) ISP, ISP1, ISM, ISM1
wwm_snl4.F90:143:            WRITE(STAT%FHNDL,*) ISLOW, ISHGH, ISCLW, ISCHG, IDLOW, IDHGH
wwm_snl4.F90:144:            WRITE(STAT%FHNDL,*) XIS
wwm_snl4.F90:145:            WRITE(STAT%FHNDL,*) AWG1, AWG2, AWG3, AWG4
wwm_snl4.F90:146:            WRITE(STAT%FHNDL,*) AWG5, AWG6, AWG7, AWG8
wwm_snl4.F90:147:            WRITE(STAT%FHNDL,*) '---------------------------------------'
wwm_snl4_tsa.F90:902:         write(*,'('' error '',i10,'' from gridset; bail'')') ierr_gr
wwm_snl4_tsa.F90:919:!wrt    write(103,903) 'kref2(npts,NTH,NZZ)'
wwm_snl4_tsa.F90:922:!wrt       write(103,904) ipt, iang
wwm_snl4_tsa.F90:923:!wrt       write(103,940) (kref2(ipt,iang,iz), iz=1,NZZ)
wwm_snl4_tsa.F90:931:!wrt    write(104,903) 'kref4(npts,NTH,NZZ)'
wwm_snl4_tsa.F90:934:!wrt       write(104,904) ipt, iang
wwm_snl4_tsa.F90:935:!wrt       write(104,940) (kref4(ipt,iang,iz), iz=1,NZZ)
wwm_snl4_tsa.F90:945:!wrt    write(105,903) 'jref2(npts,NTH,NZZ)'
wwm_snl4_tsa.F90:948:!wrt       write(105,904) ipt, iang
wwm_snl4_tsa.F90:949:!wrt       write(105,940) (jref2(ipt,iang,iz), iz=1,NZZ)
wwm_snl4_tsa.F90:957:!wrt    write(106,903) 'jref4(npts,NTH,NZZ)'
wwm_snl4_tsa.F90:960:!wrt       write(106,904) ipt, iang
wwm_snl4_tsa.F90:961:!wrt       write(106,940) (jref4(ipt,iang,iz), iz=1,NZZ)
wwm_snl4_tsa.F90:971:!wrt    write(107,903) 'wtk2(npts,NTH,NZZ)'
wwm_snl4_tsa.F90:974:!wrt       write(107,904) ipt, iang
wwm_snl4_tsa.F90:975:!wrt       write(107,950) (wtk2(ipt,iang,iz), iz=1,NZZ)
wwm_snl4_tsa.F90:983:!wrt    write(108,903) 'wtk4(npts,NTH,NZZ)'
wwm_snl4_tsa.F90:986:!wrt       write(108,904) ipt, iang
wwm_snl4_tsa.F90:987:!wrt       write(108,950) (wtk4(ipt,iang,iz), iz=1,NZZ)
wwm_snl4_tsa.F90:997:!wrt    write(109,903) 'wta2(npts,NTH,NZZ)'
wwm_snl4_tsa.F90:1000:!wrt       write(109,904) ipt, iang
wwm_snl4_tsa.F90:1001:!wrt       write(109,950) (wta2(ipt,iang,iz), iz=1,NZZ)
wwm_snl4_tsa.F90:1009:!wrt    write(110,903) 'wta4(npts,NTH,NZZ)'
wwm_snl4_tsa.F90:1012:!wrt       write(110,904) ipt, iang
wwm_snl4_tsa.F90:1013:!wrt       write(110,950) (wta4(ipt,iang,iz), iz=1,NZZ)
wwm_snl4_tsa.F90:1023:!wrt    write(111,903) 'tfac2(npts,NTH,NZZ)'
wwm_snl4_tsa.F90:1026:!wrt       write(111,904) ipt, iang
wwm_snl4_tsa.F90:1027:!wrt       write(111,950) (tfac2(ipt,iang,iz), iz=1,NZZ)
wwm_snl4_tsa.F90:1035:!wrt    write(112,903) 'tfac4(npts,NTH,NZZ)'
wwm_snl4_tsa.F90:1038:!wrt       write(112,904) ipt, iang
wwm_snl4_tsa.F90:1039:!wrt       write(112,950) (tfac4(ipt,iang,iz), iz=1,NZZ)
wwm_snl4_tsa.F90:1049:!wrt    write(113,903) 'grad(npts,NTH,NZZ)'
wwm_snl4_tsa.F90:1052:!wrt       write(113,904) ipt, iang
wwm_snl4_tsa.F90:1053:!wrt       write(113,950) (grad(ipt,iang,iz), iz=1,NZZ)
wwm_snl4_tsa.F90:1072:!b      write(123,903) 'ef2(NK,NTH) + 2 peaks'
wwm_snl4_tsa.F90:1077:!b      write(124,903) 'dens(NK,NTH) + 2 peaks'
wwm_snl4_tsa.F90:1261:        write(6,206) nbins, npk,fpk,e1max,e1sum                       !*  <<<<<
wwm_snl4_tsa.F90:1349:          write(6,216) npk2,fpk2,e1max2,nfs
wwm_snl4_tsa.F90:1352:          write(6,217) npk,fpk,e1max,h1sig
wwm_snl4_tsa.F90:1362:          write(6,218) npk,fpk,e1max,nfs
wwm_snl4_tsa.F90:1365:          write(6,219) npk2,fpk2,e1max2,h1sig
wwm_snl4_tsa.F90:1409:        write(6,207) npk,fpk,e1max,h1sig
wwm_snl4_tsa.F90:1436:        write(6,208) npk,fpk,e1max,h1sig
wwm_snl4_tsa.F90:1454:        write(6,209) npk2,fpk2,e1max2,nfs
wwm_snl4_tsa.F90:1514:      write(6,418) dens2ov1
wwm_snl4_tsa.F90:1527:!b        write(123,210) npk,fpk,e1max, npk2,fpk2,e1max2, h1sig
wwm_snl4_tsa.F90:1529:!b          write(123,211) (ef2(irng,iang), irng=1,nrng)
wwm_snl4_tsa.F90:1534:!b        write(124,210) npk,fpk,e1max, npk2,fpk2,e1max2, h1sig
wwm_snl4_tsa.F90:1536:!b          write(124,211) (dens(irng,iang), irng=1,nrng)
wwm_snl4_tsa.F90:1825:!prt  write(6,500)
wwm_snl4_tsa.F90:1826:!prt  write(6,501)
wwm_snl4_tsa.F90:1827:!prt  write(6,701) fa1max, fa1sum
wwm_snl4_tsa.F90:1828:!prt  write(6,702) fa2max, fa2sum
wwm_snl4_tsa.F90:1829:!prt  write(6,703) fa3max, fa3sum
wwm_snl4_tsa.F90:1830:!prt  write(6,704) fa4max, fa4sum
wwm_snl4_tsa.F90:1831:!prt  write(6,500)
wwm_snl4_tsa.F90:1832:!prt  write(6,705) fe1max, fe1sum
wwm_snl4_tsa.F90:1833:!prt  write(6,706) fe2max, fe2sum
wwm_snl4_tsa.F90:1834:!prt  write(6,707) fe3max, fe3sum
wwm_snl4_tsa.F90:1835:!prt  write(6,708) fe4max, fe4sum
wwm_snl4_tsa.F90:1836:!prt  write(6,511)
wwm_snl4_tsa.F90:1838:!prt  write(6,500)
wwm_snl4_tsa.F90:1839:!prt  write(6,502)
wwm_snl4_tsa.F90:1840:!prt  write(6,701) ta1max, ta1sum
wwm_snl4_tsa.F90:1841:!prt  write(6,702) ta2max, ta2sum
wwm_snl4_tsa.F90:1842:!prt  write(6,703) ta3max, ta3sum
wwm_snl4_tsa.F90:1843:!prt  write(6,704) ta4max, ta4sum
wwm_snl4_tsa.F90:1844:!prt  write(6,500)
wwm_snl4_tsa.F90:1845:!prt  write(6,705) te1max, te1sum
wwm_snl4_tsa.F90:1846:!prt  write(6,706) te2max, te2sum
wwm_snl4_tsa.F90:1847:!prt  write(6,707) te3max, te3sum
wwm_snl4_tsa.F90:1848:!prt  write(6,708) te4max, te4sum
wwm_snl4_tsa.F90:1849:!prt  write(6,511)
wwm_snl4_tsa.F90:1879:!a    write(80,'(''ZONE,T="Energy density",I=35,J=36,F=point'')')
wwm_snl4_tsa.F90:1882:!a       write(80,179)          frqa(ifr),angl(iang)*deg,ef2(ifr,iang)
wwm_snl4_tsa.F90:1890:!a    write(11,111)
wwm_snl4_tsa.F90:1895:!a      write(11,112) ifr, frqa(ifr)
wwm_snl4_tsa.F90:1903:!a    write(12,121)
wwm_snl4_tsa.F90:1906:!a      write(12,122) iang, angl(iang)*deg
wwm_snl4_tsa.F90:1909:!a    write(12,122) 37, 360.0
wwm_snl4_tsa.F90:1916:!a    write(40,*)'% Cartesian 2D Energy Density in (f,theta), I=35,J=37'
wwm_snl4_tsa.F90:1917:!b    write(40,*)'% Cartesian 2D Energy Density in (f,theta), I=28,J=37'
wwm_snl4_tsa.F90:1918:!a    write(40,*)'% 2D Enenrgy max = ', e2max
wwm_snl4_tsa.F90:1920:!a      write(40,182) ( ef2(ifr,iang), ifr=1,nrng )
wwm_snl4_tsa.F90:1921:!b      write(40,182) ( ef2(ifr,iang), ifr=1,28 )
wwm_snl4_tsa.F90:1924:!a    write(40,182) ( ef2(ifr,1), ifr=1,nrng )
wwm_snl4_tsa.F90:1925:!b    write(40,182) ( ef2(ifr,1), ifr=1,28 )
wwm_snl4_tsa.F90:1931:!a    write(81,'(''ZONE,T="Snl for energy",I=35,J=36,F=POINT'')')
wwm_snl4_tsa.F90:1934:!a       write(81,179)                                             !* &
wwm_snl4_tsa.F90:1944:!a    write(41,*)'% Polar FBI 2D Snl in (k,theta) Energy Unit,I=35,J=37'
wwm_snl4_tsa.F90:1945:!b    write(41,*)'% Polar FBI 2D Snl in (k,theta) Energy Unit,I=28,J=37'
wwm_snl4_tsa.F90:1946:!a    write(41,*)'% FBI 2D Snl max = ', fe1max
wwm_snl4_tsa.F90:1948:!a      write(41,182) ( fbie1(ifr,iang), ifr=1,nrng )
wwm_snl4_tsa.F90:1949:!b      write(41,182) ( fbie1(ifr,iang), ifr=1,28 )
wwm_snl4_tsa.F90:1951:!a    write(41,182) ( fbie1(ifr,1), ifr=1,nrng )
wwm_snl4_tsa.F90:1952:!b    write(41,182) ( fbie1(ifr,1), ifr=1,28 )
wwm_snl4_tsa.F90:1958:!a    write(84,'(''ZONE,T="TSA for energy",I=35,J=36,F=POINT'')')
wwm_snl4_tsa.F90:1961:!a       write(84,179) frqa(ifr),angl(iang)*deg, twopi*frqa(ifr)*  !* &
wwm_snl4_tsa.F90:1970:!a    write(44,*)'% Polar TSA 2D Snl in (k,theta) Energy Unit,I=35,J=37'
wwm_snl4_tsa.F90:1971:!b    write(44,*)'% Polar TSA 2D Snl in (k,theta) Energy Unit,I=28,J=37'
wwm_snl4_tsa.F90:1972:!a    write(44,*)'% TSA 2D Snl max = ', te1max
wwm_snl4_tsa.F90:1974:!a      write(44,182) ( tsae1(ifr,iang), ifr=1,nrng )
wwm_snl4_tsa.F90:1975:!b      write(44,182) ( tsae1(ifr,iang), ifr=1,28 )
wwm_snl4_tsa.F90:1977:!a    write(44,182) ( tsae1(ifr,1), ifr=1,nrng )
wwm_snl4_tsa.F90:1978:!b    write(44,182) ( tsae1(ifr,1), ifr=1,28 )
wwm_snl4_tsa.F90:1985:!a    write(85,'(''ZONE,T="Energy density-norm",I=35,J=36,F=point'')')
wwm_snl4_tsa.F90:1988:!a       write(85,179) frqa(ifr)/frqa(14),angl(iang)*deg,ef2(ifr,iang)
wwm_snl4_tsa.F90:2305:!wrt    write(125,902) 'som1,som2,som3 from cplshr at irng,krng,izz,kang,ipt'
wwm_snl4_tsa.F90:2308:!wrt    write(126,902) 'domsq23 from cplshr at irng,krng,izz,kang,ipt'
wwm_snl4_tsa.F90:2311:!wrt    write(127,902) 'sumom from cplshr at irng,krng,izz,kang,ipt'
wwm_snl4_tsa.F90:3562:!wrt    write(125,909) irng,krng,izz,kang,ipt,                        &
wwm_snl4_tsa.F90:3592:!wrt    write(126,919) irng,krng,izz,kang,ipt, sngl(domsq23)
wwm_snl4_tsa.F90:3636:!wrt  write(127,929) irng,krng,izz,kang,ipt, sngl(sumom)
wwm_snl4_tsa.F90:3917:!b      write(6,105) m, q(m)
wwm_snl4_tsa.F90:4721:!b      write(1,901)  (mat1(irng,iang), irng=1,nrng)
wwm_snl4_tsa.F90:4722:!b      write(2,901)  (mat2(irng,iang), irng=1,nrng)
wwm_snonlin.F90:139:!      WRITE(111117,'(I10,10F15.8)') IG, SUM(F), SUM(FL),SUM(SL), AKMEAN
wwm_snonlin.F90:141:!        WRITE(111117,'(5I10,10F15.8)') ISHALLO, ISNONLIN, MLSTHG, & 
wwm_snonlin.F90:144:!      WRITE(111117,'(5I10,10F15.8)') SUM(INLCOEF), SUM(IKP), SUM(IKP1), & 
wwm_snonlin.F90:146:!      WRITE(111117,'(4I10,10F20.8)') SUM(K1W), SUM(K2W), SUM(K11W), &
wwm_snonlin.F90:331:!              WRITE(111117,'(3I10,5F20.10)')KH,K,MC,FTAIL,RNLCOEF(1,MC)
wwm_snonlin.F90:332:!              WRITE(111117,'(6F20.10)') FAD1,FAD2,FCEN,DELAD(IJ)
wwm_snonlin.F90:333:!              WRITE(111117,'(5F20.15)') DELAM(IJ),FIJ
wwm_snonlin.F90:334:!              WRITE(111117,'(5F30.20)') FCEN, FTEMP(IJ), FIJ, FTAIL
wwm_snonlin.F90:528:!            WRITE(111117,'(I10,4F30.25)') MC, &
wwm_snonlin.F90:550:!        WRITE(111117,'(2F30.25)') & 
wwm_snonlin.F90:552:!        WRITE(111117,*) 'NOW THE FULL THING'
wwm_snonlin.F90:555:!            WRITE(111117,'(2I10,2F30.25)') &
wwm_sourceterms.F90:167:                 WRITE(DBG%FHNDL,*) 'NO ST42 or ST41 chosen but MESIN == 1'
wwm_sourceterms.F90:180:                   WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT NORMAL', IP, '   DUE TO SIN', SUM(IMATRA), SUM(IMATDA)
wwm_sourceterms.F90:217:             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   DUE TO SIN', SUM(IMATRA), SUM(IMATDA)
wwm_sourceterms.F90:274:             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   DUE TO SNL4'
wwm_sourceterms.F90:384:             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   DUE TO SBF' 
wwm_sourceterms.F90:418:             WRITE(DBG%FHNDL,*) 'NO ST42 or ST41 chosen but MESIN == 1'
wwm_sourceterms.F90:435:!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----SOURCE TIMINGS-----'
wwm_sourceterms.F90:436:!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPARATION        ', TIME2-TIME1
wwm_sourceterms.F90:437:!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SIN                ', TIME3-TIME2
wwm_sourceterms.F90:438:!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SDS                ', TIME4-TIME3
wwm_sourceterms.F90:439:!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SNL4               ', TIME5-TIME4
wwm_sourceterms.F90:440:!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SNL3               ', TIME6-TIME5
wwm_sourceterms.F90:441:!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SBR                ', TIME7-TIME6
wwm_sourceterms.F90:442:!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SBF                ', TIME8-TIME7
wwm_sourceterms.F90:443:!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'RECALC             ', TIME9-TIME8
wwm_sourceterms.F90:444:!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'TOTAL              ', TIME9-TIME1
wwm_sourceterms.F90:445:!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '------END-TIMINGS-  ---'
wwm_sourceterms.F90:559:!        WRITE(300,*) 'IP=', IP, ' HMAX=', HMAX(IP), ' DEP=', DEP(IP)
wwm_sourceterms.F90:560:!        WRITE(300,*) '   ', IP, ' EMAX=', EMAX, ' ETOT=', ETOT
wwm_sourceterms.F90:561:!        WRITE(300,*) '   ', IP, ' HS=', HS, ' BRHD=', BRHD
wwm_sourceterms.F90:564:          WRITE(300,*) '   break XP=', XP(IP)
wwm_sparskit.F90:2858:!                   write(DBG%FHNDL,*) ii, jw, jj
wwm_sparskit.F90:2994:!                    write(DBG%FHNDL,*) ii, jw, jj
wwm_sparskit.F90:3000:              write(DBG%FHNDL,*) 'zero pivot'
wwm_specint.F90:72:               WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   IN SOURCE TERM INTEGRATION'
wwm_specint.F90:173:                 WRITE(111112,'(A10,I10)') 'BEFORE', IP
wwm_specint.F90:174:                 WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(IP,:,:))
wwm_specint.F90:175:                 WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(IP,:,:))
wwm_specint.F90:176:                 WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
wwm_specint.F90:177:                 WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
wwm_specint.F90:178:                 WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
wwm_specint.F90:179:                 WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
wwm_specint.F90:180:                 WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
wwm_specint.F90:181:                 WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
wwm_specint.F90:182:                 WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
wwm_specint.F90:183:                 WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
wwm_specint.F90:184:                 WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
wwm_specint.F90:185:                 WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
wwm_specint.F90:186:                 WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
wwm_specint.F90:187:                 WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(IP,:,:))
wwm_specint.F90:188:                 WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(IP,:))
wwm_specint.F90:198:                 WRITE(111112,'(A10,I10)') 'AFTER', IP
wwm_specint.F90:199:                 WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(IP,:,:))
wwm_specint.F90:200:                 WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(IP,:,:))
wwm_specint.F90:201:                 WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
wwm_specint.F90:202:                 WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
wwm_specint.F90:203:                 WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
wwm_specint.F90:204:                 WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
wwm_specint.F90:205:                 WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
wwm_specint.F90:206:                 WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
wwm_specint.F90:207:                 WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
wwm_specint.F90:208:                 WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
wwm_specint.F90:209:                 WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
wwm_specint.F90:210:                 WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
wwm_specint.F90:211:                 WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
wwm_specint.F90:212:                 WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(IP,:,:))
wwm_specint.F90:213:                 WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(IP,:))
wwm_specint.F90:231:                   IF (NEWDAC/MAXDAC .gt. one) WRITE(*,*) ONE/MAX(ONE,NEWDAC/MAXDAC), NEWDAC/MAXDAC
wwm_specint.F90:278:                   WRITE(111112,'(A10,I10)') 'BEFORE', IP
wwm_specint.F90:279:                   WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(IP,:,:))
wwm_specint.F90:280:                   WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(IP,:,:))
wwm_specint.F90:281:                   WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
wwm_specint.F90:282:                   WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
wwm_specint.F90:283:                   WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
wwm_specint.F90:284:                   WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
wwm_specint.F90:285:                   WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
wwm_specint.F90:286:                   WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
wwm_specint.F90:287:                   WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
wwm_specint.F90:288:                   WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
wwm_specint.F90:289:                   WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
wwm_specint.F90:290:                   WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
wwm_specint.F90:291:                   WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
wwm_specint.F90:292:                   WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(IP,:,:))
wwm_specint.F90:293:                   WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(1,:))
wwm_specint.F90:303:                   WRITE(111112,'(A10,I10)') 'AFTER', IP
wwm_specint.F90:304:                   WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(IP,:,:))
wwm_specint.F90:305:                   WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(IP,:,:))
wwm_specint.F90:306:                   WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
wwm_specint.F90:307:                   WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
wwm_specint.F90:308:                   WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
wwm_specint.F90:309:                   WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
wwm_specint.F90:310:                   WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
wwm_specint.F90:311:                   WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
wwm_specint.F90:312:                   WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
wwm_specint.F90:313:                   WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
wwm_specint.F90:314:                   WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
wwm_specint.F90:315:                   WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
wwm_specint.F90:316:                   WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
wwm_specint.F90:317:                   WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(IP,:,:))
wwm_specint.F90:318:                   WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(IP,:))
wwm_specint.F90:333:         !WRITE(*,*) SUM(IMATRAA), SUM(IMATDAA)
wwm_specint.F90:366:               IF (IP == TESTNODE) WRITE(*,'(A20,3F15.8)') 'POST BEFORE', SUM(SL(IP,:,:)), SUM(FL3(IP,:,:)),  SUM(FL(IP,:,:))
wwm_specint.F90:374:               IF (IP == TESTNODE) WRITE(*,'(A20,3F15.8)') 'POST AFTER', SUM(SL(IP,:,:)), SUM(FL3(IP,:,:)),  SUM(FL(IP,:,:))
wwm_specint.F90:470:!         write(*,'(A10,2I10,L10,I10,2F15.6)') 'after', ip, iobp(ip), limiter, iselect, sum(acloc), sum(imatra)
wwm_specint.F90:520:         !WRITE(*,*) '1 RK-TVD', SUM(ACOLD), SUM(ACLOC)
wwm_specint.F90:533:         !WRITE(*,*) '2 RK-TVD', SUM(ACOLD), SUM(ACLOC)
wwm_specint.F90:546:         !WRITE(*,*) '3 RK-TVD', SUM(ACOLD), SUM(ACLOC)
wwm_specparam.F90:208:!         WRITE(*,*) FP, ETOTF3, ETOTF4
wwm_specparam.F90:1019:              WRITE(*,*) ETOT_SKD, tmp(is), tmp(is-1), ACLOC(IS,ID), WK(IP,IS)*DEP(IP), CALLFROM
wwm_specparam.F90:1026:              WRITE(*,*) ETOT_SKD, tmp(msc),CALLFROM
wwm_specparam.F90:1131:         !WRITE(*,'(9F15.4)') DEPLOC,SUM(ACLOC),CURTXYLOC,SUM(WKLOC), ETOT_SKD, ETOT_SKDSIG 
wwm_specparam.F90:1390:       !WRITE(*,'(11F15.4)') FPP, KPP, CGPP, WKDEPP, WNPP, CPP, TPP, LPP, PEAKDM, PEAKFF, PEAKDSPR 
wwm_stress.F90:100:!        WRITE(STAT%FHNDL,*) 'ALLOCATED STRESS TABLE', ITAUMAX, JUMAx, JPLEVT
wwm_stress.F90:109:      IF (LOUTWAM) WRITE(111111,'(A20)') 'STRESS'
wwm_stress.F90:115:      IF (LOUTWAM) WRITE(111111,'(A30,I10)') 'TOTAL NUMBER OF ENTRIES -- STRESS --', &
wwm_stress.F90:118:      WRITE(5011) DELU, DELTAUW
wwm_stress.F90:128:          WRITE(IU06,*)' '
wwm_stress.F90:129:          WRITE(IU06,*)' STRESS FOR LEVEL HEIGHT ',XL,' m'
wwm_stress.F90:132:        IF (LOUTWAM) WRITE(111111,'(3I10,F15.8)') JL, JPLEVT, JPLEVC, XL
wwm_stress.F90:155:!            WRITE(111111,'(10F15.8)') I,J,JL,TAUT(I,J,JL),UTOP,USTOLD,TAUOLD
wwm_stress.F90:171:      WRITE(5011) TAUT
wwm_stress.F90:172:      IF (LOUTWAM) WRITE(111111,'(3F15.6)') DELU,DELTAUW,SUM(TAUT)
wwm_stresso.F90:133:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,*) '------- STARTING STRESSO ---------'
wwm_stresso.F90:155:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(10F15.8)') CONST, ROG, ROWATER, G, CONST0(IJ), SUM(ZPIROFR)
wwm_stresso.F90:166:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(10F15.8)') SCDFM, CONSTFM(IJ,MIJ(IJ)) 
wwm_stresso.F90:208:        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(2F15.8)') XSTRESS(IJ), YSTRESS(IJ)
wwm_stresso.F90:264:          IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(4F15.8)') USDIRP(IJ), COSW, TEMP1(IJ), TEMP2(IJ)
wwm_stresso.F90:292:          IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(A10,2F20.10,2I10)') 'T1', XI, XJ, I, J
wwm_stresso.F90:293:          IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(A10,4F20.10)') 'T2', DELI2, DELI1, DELJ2, DELJ1
wwm_stresso.F90:294:          IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(A10,4F20.10)') 'T3', CONST0(IJ), TEMP1(IJ), UST2(IJ), TAU1
wwm_stresso.F90:340:         IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,'(4F15.8)') XSTRESS(IJ), YSTRESS(IJ) , TAUW(IJ), TAUHF(IJ)
wwm_stresso.F90:356:      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111116,*) '------- END OF STRESSO ---------'
wwm_tauhf.F90:91:      IF (LOUTWAM) WRITE(111111,'(A20)') 'TAUHF'
wwm_tauhf.F90:93:      IF (LOUTWAM) WRITE(111111,'(10F15.10)') ALPHAM, DELUST, DELALP, DELTAIL, CONST1, SUM(W)
wwm_tauhf.F90:101:      IF (LOUTWAM) WRITE(111111,'(A20,I10)') 'TOTAL NUMBER OF ENTRIES', ML*IALPHA*IUSTAR
wwm_tauhf.F90:138:!            WRITE(111111,'(3I10,10F15.7)') L,K,J,DELY,ZMU,TAUHFT(K,L,M)
wwm_tauhf.F90:195:      IF (LOUTWAM) WRITE(111111,'(A10,I10)') 'IPHYS=', IPHYS
wwm_tauhf.F90:197:        WRITE(5011) DELALP, DELUST, DELTAIL
wwm_tauhf.F90:198:        WRITE(5011) TAUHFT
wwm_tauhf.F90:199:        IF (LOUTWAM) WRITE(111111,'(F20.10)') SUM(TAUHFT)
wwm_tauhf.F90:201:        WRITE(5011) DELALP, DELUST, DELTAIL
wwm_tauhf.F90:202:        WRITE(5011) TAUHFT, TAUHFT2, TAUW
wwm_tauhf.F90:203:        IF (LOUTWAM) WRITE(111111,'(3F20.10)') DELTAIL, SUM(TAUHFT), SUM(TAUHFT2) 
wwm_tauhf.F90:207:!        WRITE(100000,*) M, FR(M)
wwm_tauhf.F90:213:!            WRITE(100002,*) K,L,M,TAUHFT(K,L,M)
wwm_wave_setup.F90:39:!        WRITE(700,*) 'ISS=', ISS, 'INCR=', DS_INCR(ISS), 'diff=', SPSIG(ISS) - SPSIG(ISS-1)
wwm_wave_setup.F90:64:!        WRITE(700,*) 'IP=', IP, 'HS=', eHS, 'RXX/RYY=', eRXX, eRYY
wwm_wave_setup.F90:197:!            WRITE(200+myrank,*) 'UGRAD=', UGRAD, 'VGRAD=', VGRAD
wwm_wave_setup.F90:198:!            WRITE(200+myrank,*) 'UGRAD1=', UGRAD1, 'VGRAD1=', VGRAD1
wwm_wave_setup.F90:199:!            WRITE(200+myrank,*) 'I1=', I1, ' K=', K
wwm_wave_setup.F90:200:!            WRITE(200+myrank,*) 'I1=', I1, ' IDX=', IDX, 'eScal=', eScal
wwm_wave_setup.F90:202:!            WRITE(200+myrank,*) 'IPp=', IPp, ' JPp=', JPp
wwm_wave_setup.F90:203:!            WRITE(200+myrank,*) '            -  -  -  -  -'
wwm_wave_setup.F90:209:!        WRITE(200+myrank,*) '--------------------------------------'
wwm_wave_setup.F90:239:!            WRITE(200+myrank,*) 'aspar(J1)=', ASPAR(J1)
wwm_wave_setup.F90:240:!            WRITE(200+myrank,*) 'aspar(J2)=', ASPAR(J2)
wwm_wave_setup.F90:286:            WRITE(*,*) 'IP=', IP, 'J=', J, ' nbM=', nbM
wwm_wave_setup.F90:291:      WRITE(200 + myrank,*) 'Symmetry error=', eSum
wwm_wave_setup.F90:380:      WRITE(200+myrank,*) 'sum(V_R)=', sum(V_R)
wwm_wave_setup.F90:381:      WRITE(200+myrank,*) 'sum(V_Z)=', sum(V_Z)
wwm_wave_setup.F90:382:      WRITE(200+myrank,*) 'Before loop, |B|=', eNorm
wwm_wave_setup.F90:388:        WRITE(200+myrank,*) 'nbIter=', nbIter
wwm_wave_setup.F90:389:        WRITE(200+myrank,*) 'Before call to WAVE_SETUP_APPLY_FCT'
wwm_wave_setup.F90:394:        WRITE(200+myrank,*) 'After call to WAVE_SETUP_APPLY_FCT'
wwm_wave_setup.F90:400:        WRITE(200+myrank,*) 'sum(V_P)=', sum(V_P)
wwm_wave_setup.F90:401:        WRITE(200+myrank,*) 'sum(V_Y)=', sum(V_Y)
wwm_wave_setup.F90:402:        WRITE(200+myrank,*) 'h2=', h2
wwm_wave_setup.F90:403:        WRITE(200+myrank,*) 'alphaV=', alphaV
wwm_wave_setup.F90:414:        WRITE(200+myrank,*) 'nbIter=', nbIter, 'eNorm=', eNorm
wwm_wave_setup.F90:431:      WRITE(STAT%FHNDL,*) 'wave_setup nbIter=', nbIter
wwm_wave_setup.F90:434:      WRITE(200+myrank,*) 'MNP=', MNP, ' TheOut:'
wwm_wave_setup.F90:436:        WRITE(200+myrank,*) 'IP=', IP, ' setup=', TheOut(IP)
wwm_wave_setup.F90:533:      WRITE(200 + myrank,*) 'WAVE_SETUP_COMPUTATION, step 1'
wwm_wave_setup.F90:539:      WRITE(200 + myrank,*) 'WAVE_SETUP_COMPUTATION, step 2'
wwm_wave_setup.F90:548:      WRITE(200 + myrank,*) 'sum(abs(ASPAR))=', sum(abs(ASPAR))
wwm_wave_setup.F90:549:      WRITE(200 + myrank,*) 'eResidual=', eResidual
wwm_wave_setup.F90:550:      WRITE(200 + myrank,*) 'eResidual2=', eResidual2
wwm_wave_setup.F90:551:      WRITE(200 + myrank,*) 'WAVE_SETUP_COMPUTATION, step 3'
wwm_wave_setup.F90:562:      WRITE(200 + myrank,*) 'Before DEBUG statement'
wwm_wave_setup.F90:564:      WRITE(200 + myrank,*) 'After DEBUG statement'
wwm_wave_setup.F90:567:      WRITE(200 + myrank,*) 'Norm(B)=', eNorm
wwm_wave_setup.F90:571:      WRITE(200 + myrank,*) 'Norm(residual)=', eNorm
wwm_wind.F90:43:          WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'HOMOGENOUS STEADY WIND FIELD IS USED' 
wwm_wind.F90:44:          WRITE(WINDBG%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT, LWDIR
wwm_wind.F90:59:          WRITE(WINDBG%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
wwm_wind.F90:60:          WRITE(WINDBG%FHNDL,'("+TRACE...",A)')  'SPATIAL VARIABLE WIND FIELD IS USED'
wwm_wind.F90:84:            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'COMPUTING CF INTERPOLATION COEFS AND LOADING WIND_TIME_MJD'
wwm_wind.F90:133:          WRITE(WINDBG%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
wwm_wind.F90:134:          WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'NONSTATIONARY WIND FIELD IS USED        '
wwm_wind.F90:139:          WRITE(WINDBG%FHNDL,*) SEWI%BEGT, SEWI%ENDT, SEWI%ISTP, SEWI%TOTL/3600.0, SEWI%DELT
wwm_wind.F90:140:          WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'SPATIAL VARIABLE WIND FIELD IS USED'
wwm_wind.F90:141:          WRITE(WINDBG%FHNDL,*) 'IWINDFORMAT=', IWINDFORMAT
wwm_wind.F90:166:            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'SPATIAL/TEMPORAL VARIABLE WIND FIELD IS USED CF NETCDF'
wwm_wind.F90:167:            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'COMPUTING CF INTERPOLATION COEFS AND LOADING WIND_TIME_MJD'
wwm_wind.F90:181:            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'SPATIAL/TEMPORAL VARIABLE WIND FIELD IS USED CF NETCDF'
wwm_wind.F90:182:            WRITE(WINDBG%FHNDL,'("+TRACE...",A)') 'COMPUTING CF INTERPOLATION COEFS AND LOADING WIND_TIME_MJD'
wwm_wind.F90:213:      write(WINDBG%FHNDL,'("+TRACE... Done with CF init, Uwind ",F7.2,2x,F7.2)')minval(WINDXY(:,1)),maxval(WINDXY(:,1))
wwm_wind.F90:214:      write(WINDBG%FHNDL,'("+TRACE... Done with CF init, Vwind ",F7.2,2x,F7.2)')minval(WINDXY(:,2)),maxval(WINDXY(:,2))
wwm_wind.F90:236:      WRITE(WINDBG%FHNDL,*) MAIN%TMJD, SEWI%TMJD-1.E-8, MAIN%TMJD .ge. SEWI%TMJD-1.E-8, MAIN%TMJD .le. SEWI%EMJD, SEWI%EMJD
wwm_wind.F90:318:      write(WINDBG%FHNDL,'("max WINDXY:",2F7.2)')maxval(WINDXY(:,1)),maxval(WINDXY(:,2))
wwm_wind.F90:319:      write(WINDBG%FHNDL,'("min WINDXY:",2F7.2)')minval(WINDXY(:,1)),minval(WINDXY(:,2))
wwm_wind.F90:322:        WRITE(3333,*) SEWI%TMJD
wwm_wind.F90:323:        WRITE(3333,*) WINDXY(:,1)
wwm_wind.F90:324:        WRITE(3333,*) WINDXY(:,2)
wwm_wind.F90:361:        WRITE(WINDBG%FHNDL,*) 'NDT_WIND_ALL_FILES=', NDT_WIND_ALL_FILES
wwm_wind.F90:364:          WRITE(WINDBG%FHNDL,*) ' eIdx=', eIdx
wwm_wind.F90:365:          WRITE(WINDBG%FHNDL,*) ' eTime=', WIND_TIME_ALL_FILES(eIdx)
wwm_wind.F90:366:          WRITE(WINDBG%FHNDL,*) ' eTimeStr=', eTimeStr
wwm_wind.F90:411:      WRITE(WINDBG%FHNDL,*) 'KERNEL_INTERP_UV_WINDFD'
wwm_wind.F90:412:      WRITE(WINDBG%FHNDL,*) 'UWIND_FD, min/max=', minval(UWIND_FD), maxval(UWIND_FD)
wwm_wind.F90:413:      WRITE(WINDBG%FHNDL,*) 'VWIND_FD, min/max=', minval(VWIND_FD), maxval(VWIND_FD)
wwm_wind.F90:414:      WRITE(WINDBG%FHNDL,*) 'UWIND_FE, min/max=', minval(outwind(:,1)), maxval(outwind(:,1))
wwm_wind.F90:415:      WRITE(WINDBG%FHNDL,*) 'VWIND_FE, min/max=', minval(outwind(:,2)), maxval(outwind(:,2))
wwm_wind.F90:416:!      WRITE(WINDBG%FHNDL,*) 'max(CF_COEFF)=', maxval(abs(CF_COEFF))
wwm_wind.F90:417:!      WRITE(WINDBG%FHNDL,*) 'cf_scale_factor=', cf_scale_factor
wwm_wind.F90:418:!      WRITE(WINDBG%FHNDL,*) 'cf_add_offset=', cf_add_offset
wwm_wind.F90:441:      WRITE(WINDBG%FHNDL,*) 'Starting node loop for calcs of coefs'
wwm_wind.F90:457:        WRITE(WINDBG%FHNDL,*) 'min(lon)=', minval(lon)
wwm_wind.F90:458:        WRITE(WINDBG%FHNDL,*) 'max(lon)=', maxval(lon)
wwm_wind.F90:459:        WRITE(WINDBG%FHNDL,*) 'min(lat)=', minval(lat)
wwm_wind.F90:460:        WRITE(WINDBG%FHNDL,*) 'max(lat)=', maxval(lat)
wwm_wind.F90:545:              WRITE(WINDBG%FHNDL,*) 'aShift=', aShift
wwm_wind.F90:546:              WRITE(WINDBG%FHNDL,*) 'outside node IP=', I
wwm_wind.F90:547:              WRITE(WINDBG%FHNDL,*) 'eX=', eX, 'eY=', eY
wwm_wind.F90:618:      WRITE(WINDBG%FHNDL,*) ' sum(StatusUse)=', sum(StatusUse)
wwm_wind.F90:619:      WRITE(WINDBG%FHNDL,*) ' done interp calcs'
wwm_wind.F90:647:      WRITE(WINDBG%FHNDL,*) 'Time error in wind for CF'
wwm_wind.F90:648:      WRITE(WINDBG%FHNDL,*) 'MAIN % TMJD=', MAIN%TMJD
wwm_wind.F90:649:      WRITE(WINDBG%FHNDL,*) 'min(wind_time_mjd)=', minval(wind_time_mjd)
wwm_wind.F90:650:      WRITE(WINDBG%FHNDL,*) 'max(wind_time_mjd)=', maxval(wind_time_mjd)
wwm_wind.F90:688:!        WRITE(WINDBG%FHNDL,*) IT, NETCDF_FILE_NAMES(IT)
wwm_wind.F90:798:!          WRITE(WINDBG%FHNDL,*) IFILE, IT, WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE)
wwm_wind.F90:970:        WRITE(3011) TIME
wwm_wind.F90:971:        WRITE(3011) (U(IP), V(IP), H(IP), IP = 1, numLons*numLats)
wwm_wind.F90:1031:        WRITE(WINDBG%FHNDL,*) IT, NETCDF_FILE_NAMES(IT)
wwm_wind.F90:1116:      WRITE(WINDBG%FHNDL,*) 'NDT_WIND_FILE=', NDT_WIND_FILE
wwm_wind.F90:1117:      WRITE(WINDBG%FHNDL,*) 'START_TIME=', START_TIME
wwm_wind.F90:1131:        WRITE(WINDBG%FHNDL,*) 'IFILE=', IFILE
wwm_wind.F90:1140:      WRITE(WINDBG%FHNDL,*) 'NDT_WIND_ALL_FILES=', NDT_WIND_ALL_FILES
wwm_wind.F90:1154:        WRITE(WINDBG%FHNDL,*) 'IFILE=', IFILE
wwm_wind.F90:1165:              WRITE(WINDBG%FHNDL,110) IT, idx, WIND_TIME_NETCDF(IT) * 3600. * SEC2DAY, eNewTime, eTimeStr
wwm_wind.F90:1175:            WRITE(WINDBG%FHNDL,110) IT, idx, WIND_TIME_NETCDF(IT) * 3600. * SEC2DAY, eNewTime, eTimeStr
wwm_wind.F90:1182:      WRITE(WINDBG%FHNDL,*) 'WIND TIME STEP', SEWI%DELT, DAY2SEC, WIND_TIME_ALL_FILES(2), WIND_TIME_ALL_FILES(1)
wwm_wind.F90:1228:      WRITE(WINDBG%FHNDL,*) 'Begin INIT_NETCDF_NARR'
wwm_wind.F90:1229:      WRITE(WINDBG%FHNDL,*) 'MULTIPLE_IN_WIND=', MULTIPLE_IN_WIND
wwm_wind.F90:1253:      WRITE(WINDBG%FHNDL,*) 'NUM_NETCDF_FILES=', NUM_NETCDF_FILES
wwm_wind.F90:1257:        WRITE(WINDBG%FHNDL,*) 'IT=', IT, 'file=', NETCDF_FILE_NAMES(IT)
wwm_wind.F90:1276:      WRITE(WINDBG%FHNDL,*) 'NDT_WIND_FILE=', NDT_WIND_FILE
wwm_wind.F90:1324:      WRITE(WINDBG%FHNDL,*) 'OFFSET_X_WIND=', OFFSET_X_WIND
wwm_wind.F90:1325:      WRITE(WINDBG%FHNDL,*) 'OFFSET_Y_WIND=', OFFSET_Y_WIND
wwm_wind.F90:1339:      WRITE(WINDBG%FHNDL,*) 'NDT_WIND_ALL_FILES=', NDT_WIND_ALL_FILES
wwm_wind.F90:1377:        !WRITE(WINDBG%FHNDL,*) CHRDATE, START_TIME
wwm_wind.F90:1382:        WRITE(WINDBG%FHNDL,*) 'Just leaving'
wwm_wind.F90:1450:          WRITE(WINDBG%FHNDL,*) 'POINT OF THE MESH IS OUT OF THE WIND FIELD', IP, XP(IP), YP(IP)
wwm_wind.F90:1458:!         WRITE(WINDBG%FHNDL,*) 'IP=', MNP, ' sumWi=', sum(Wi)
wwm_wind.F90:1459:!         WRITE(WINDBG%FHNDL,*) 'IP=', MNP, ' minW=', minval(Wi), ' maxW=', maxval(Wi)
wwm_wind.F90:1462:      WRITE(WINDBG%FHNDL,*) 'MNP_WIND=', MNP_WIND, ' nbFail=', nbFail
wwm_wind.F90:1665:        WRITE(WINDBG%FHNDL,*) 'READ_NETCDF_CRFS IFILE=', IFILE, ' IT=', IT
wwm_wind.F90:1785:          WRITE(3011) TIME
wwm_wind.F90:1786:          WRITE(3011) (U(IP), V(IP), SQRT(U(IP)**2.+V(IP)**2.), IP = 1, NDX_WIND*NDY_WIND)
wwm_wind.F90:1868:        !WRITE(WINDBG%FHNDL,*) scale_factor
wwm_wind.F90:1869:        !WRITE(WINDBG%FHNDL,*) numLons, numLats, numTime
wwm_wind.F90:1870:        !WRITE(WINDBG%FHNDL,*) NDX_WIND, NDY_WIND
wwm_wind.F90:1966:          WRITE(3011) TIME
wwm_wind.F90:1967:          WRITE(3011) (UWND_NARR(IP), VWND_NARR(IP), SQRT(UWND_NARR(IP)**2.+VWND_NARR(IP)**2.), IP = 1, NDX_WIND*NDY_WIND)
wwm_wind.F90:1994:        WRITE(WINDBG%FHNDL,*) 'ErrorCoord=', ErrorCoord
wwm_wind.F90:2023:        WRITE(WINDBG%FHNDL,*) 'END OF RUN'
wwm_wind.F90:2024:        WRITE(WINDBG%FHNDL,*) 'WIND START TIME is outside CF wind_time range!'
wwm_wind.F90:2026:        WRITE(WINDBG%FHNDL,*) 'SEWI%BMJD=', SEWI%BMJD, ' date=', eTimeStr
wwm_wind.F90:2028:        WRITE(WINDBG%FHNDL,*) 'SEWI%EMJD=', SEWI%EMJD, ' date=', eTimeStr
wwm_wind.F90:2030:        WRITE(WINDBG%FHNDL,*) 'min(WIND_TIME_MJD)=', minval(WIND_TIME_MJD), ' date=', eTimeStr
wwm_wind.F90:2032:        WRITE(WINDBG%FHNDL,*) 'max(WIND_TIME_MJD)=', maxval(WIND_TIME_MJD), ' date=', eTimeStr
wwm_wind.F90:2034:        WRITE(wwmerr, *) 'Error in WIND_TIME_MJD 1, read ', TRIM(WINDBG%FNAME)
wwm_wind.F90:2038:        WRITE(WINDBG%FHNDL,*) 'END OF RUN'
wwm_wind.F90:2039:        WRITE(WINDBG%FHNDL,*) 'WIND END TIME is outside CF wind_time range!'
wwm_wind.F90:2041:        WRITE(WINDBG%FHNDL,*) 'SEWI%BMJD=', SEWI%BMJD, ' date=', eTimeStr
wwm_wind.F90:2043:        WRITE(WINDBG%FHNDL,*) 'SEWI%EMJD=', SEWI%EMJD, ' date=', eTimeStr
wwm_wind.F90:2045:        WRITE(WINDBG%FHNDL,*) 'min(WIND_TIME_MJD)=', minval(WIND_TIME_MJD), ' date=', eTimeStr
wwm_wind.F90:2047:        WRITE(WINDBG%FHNDL,*) 'max(WIND_TIME_MJD)=', maxval(WIND_TIME_MJD), ' date=', eTimeStr
wwm_wind.F90:2049:        WRITE(wwmerr, *) 'Error in WIND_TIME_MJD 2, read ', TRIM(WINDBG%FNAME)
wwm_wind.F90:2141:      WRITE(WINDBG%FHNDL,*) 'WindTimeStr=', TRIM(WindTimeStr)
wwm_wind.F90:2142:      WRITE(WINDBG%FHNDL,*) 'Checking for scale_factor'
wwm_wind.F90:2148:        WRITE(WINDBG%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
wwm_wind.F90:2151:      WRITE(WINDBG%FHNDL,*) 'cf_scale_factor=', cf_scale_factor
wwm_wind.F90:2154:      WRITE(WINDBG%FHNDL,*) 'Checking for add_offset'
wwm_wind.F90:2159:        WRITE(WINDBG%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
wwm_wind.F90:2162:      WRITE(WINDBG%FHNDL,*) 'cf_add_offset=', cf_add_offset
wwm_wind.F90:2171:      WRITE(WINDBG%FHNDL,*) 'Xname=', TRIM(Xname)
wwm_wind.F90:2172:      WRITE(WINDBG%FHNDL,*) 'Yname=', TRIM(Yname)
wwm_wind.F90:2189:      WRITE(WINDBG%FHNDL,*) 'NDX_WIND_FD=', NDX_WIND_FD
wwm_wind.F90:2190:      WRITE(WINDBG%FHNDL,*) 'NYX_WIND_FD=', NDY_WIND_FD
wwm_wind.F90:2216:      WRITE(WINDBG%FHNDL,*) 'eStrUnitTime=', TRIM(eStrUnitTime)
wwm_wind.F90:2217:      WRITE(WINDBG%FHNDL,*) 'eTimeStart=', eTimeStart
wwm_wind.F90:2313:      WRITE(WINDBG%FHNDL,*) 'READ_DIRECT_NETCDF_CF'
wwm_wind.F90:2314:      WRITE(WINDBG%FHNDL,*) 'RECORD_IN=', RECORD_IN
wwm_wind.F90:2315:      WRITE(WINDBG%FHNDL,*) 'UWIND_FD, min/max=', minval(UWIND_FD), maxval(UWIND_FD)
wwm_wind.F90:2316:      WRITE(WINDBG%FHNDL,*) 'VWIND_FD, min/max=', minval(VWIND_FD), maxval(VWIND_FD)
wwm_wind.F90:2317:      WRITE(WINDBG%FHNDL,*) 'UWIND_FE, min/max=', minval(outwind(:,1)), maxval(outwind(:,1))
wwm_wind.F90:2318:      WRITE(WINDBG%FHNDL,*) 'VWIND_FE, min/max=', minval(outwind(:,2)), maxval(outwind(:,2))
wwm_wind.F90:2350:      WRITE(WINDBG%FHNDL,*) 'variable used for time=', TRIM(WindTimeStr)
wwm_wind.F90:2356:        WRITE(WINDBG%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
wwm_wind.F90:2359:      WRITE(WINDBG%FHNDL,*) 'cf_scale_factor=', cf_scale_factor
wwm_wind.F90:2365:        WRITE(WINDBG%FHNDL,*) 'CHRERR=', TRIM(CHRERR)
wwm_wind.F90:2368:      WRITE(WINDBG%FHNDL,*) 'cf_add_offset=', cf_add_offset
wwm_wind.F90:2382:      WRITE(WINDBG%FHNDL,*) 'eTimeStart=', eTimeStart
wwm_wind.F90:2436:      WRITE(WINDBG%FHNDL,*) 'NUM_GRIB_FILES=', NUM_GRIB_FILES
wwm_wind.F90:2442:        WRITE(WINDBG%FHNDL,*) IT, GRIB_FILE_NAMES(IT)
wwm_wind.F90:2451:        WRITE(WINDBG%FHNDL, *) '---------------------------------------'
wwm_wind.F90:2452:        WRITE(WINDBG%FHNDL, *) 'IT=', IT, 'file = ',  GRIB_FILE_NAMES(IT)
wwm_wind.F90:2453:        WRITE(WINDBG%FHNDL, *) 'SHIFT_WIND_TIME=', SHIFT_WIND_TIME
wwm_wind.F90:2460:        WRITE(WINDBG%FHNDL, *) 'dataDate=', dataDate
wwm_wind.F90:2471:        WRITE(WINDBG%FHNDL, *) 'stepRange=', stepRange
wwm_wind.F90:2474:          WRITE(WINDBG%FHNDL, *) 'dataTime=', dataTime
wwm_wind.F90:2483:        WRITE(WINDBG%FHNDL, *) 'IT=', IT, 'Year/m/d=', eYear, eMonth, eDay
wwm_wind.F90:2484:        WRITE(WINDBG%FHNDL, *) 'IT=', IT, 'Hour/m/s=', eHour, eMin, eSec
wwm_wind.F90:2485:        WRITE(eStrTime,10) eYear, eMonth, eDay, eHour, eMin, eSec
wwm_wind.F90:2489:        WRITE(WINDBG%FHNDL, *) 'eTimeMjd=', eTimeMjd
wwm_wind.F90:2516:          WRITE(WINDBG%FHNDL, *) 'NDX_WIND_FD=', NDX_WIND_FD
wwm_wind.F90:2517:          WRITE(WINDBG%FHNDL, *) 'NDY_WIND_FD=', NDY_WIND_FD
wwm_wind.F90:2530:          WRITE(WINDBG%FHNDL, *) 'LONGITUDE'
wwm_wind.F90:2531:          WRITE(WINDBG%FHNDL, *) 'longitudeOfFirstGridPointInDegrees=', longitudeOfFirstPointInDegrees
wwm_wind.F90:2532:          WRITE(WINDBG%FHNDL, *) 'longitudeOfLastGridPointInDegrees=', longitudeOfLastPointInDegrees
wwm_wind.F90:2533:          WRITE(WINDBG%FHNDL, *) 'LATITUDE'
wwm_wind.F90:2534:          WRITE(WINDBG%FHNDL, *) 'latitudeOfFirstGridPointInDegrees=', latitudeOfFirstPointInDegrees
wwm_wind.F90:2535:          WRITE(WINDBG%FHNDL, *) 'latitudeOfLastGridPointInDegrees=', latitudeOfLastPointInDegrees
wwm_wind.F90:2537:          WRITE(WINDBG%FHNDL, *) 'iDirectionIncrement=', iDirectionIncrement
wwm_wind.F90:2538:          WRITE(WINDBG%FHNDL, *) 'jDirectionIncrement=', jDirectionIncrement
wwm_wind.F90:2581:        WRITE(WINDBG%FHNDL,*) 'IT=', IT, 'file = ',  GRIB_FILE_NAMES(IT)
wwm_wind.F90:2584:        WRITE(WINDBG%FHNDL,*) 'n=', n
wwm_windinput.F90:72:                !write(*,'(i10,3F15.6)') ip, wind10, cd(ip), ufric(ip)
wwm_windinput.F90:127:!  WRITE(BG%FHNDL,*)  EPS_D**2, HS_W, KP_W, EXP(-EPS_T**2/EPS_D**2)
wwm_windinput.F90:132:!                  WRITE(DBG%FHNDL,*) 'ITERATION  =', I
wwm_windinput.F90:133:!                  WRITE(DBG%FHNDL,*) 'wind10     =', wind10
wwm_windinput.F90:134:!                  WRITE(DBG%FHNDL,*) 'fU10       =', fU10
wwm_windinput.F90:135:!                  WRITE(DBG%FHNDL,*) 'z0_t       =', z0_t
wwm_windinput.F90:136:!                  WRITE(DBG%FHNDL,*) 'z0         =', z0(ip)
wwm_windinput.F90:137:!                  WRITE(DBG%FHNDL,*) 'Hs         =', HS_w, 'L =', Lur, 'TP =', TP_w
wwm_windinput.F90:138:!                  WRITE(DBG%FHNDL,*) 'Ulur       =' ,ULur, 'Ur =', Ur, 'EPS_D =', EPS_D
wwm_windinput.F90:139:!                  WRITE(DBG%FHNDL,*) 'T_ds & T_t =', TAUW(ip), TAUHF(ip)
wwm_windinput.F90:140:!                  WRITE(DBG%FHNDL,*) 'UFRIC      =', UFRIC1, UFRIC2
wwm_windinput.F90:154:                !WRITE(DBG%FHNDL,*) 'ITERATION  =', I
wwm_windinput.F90:155:                !WRITE(DBG%FHNDL,*) 'wind10     =', wind10
wwm_windinput.F90:156:                !WRITE(DBG%FHNDL,*) 'fU10       =', fU10
wwm_windinput.F90:157:                !WRITE(DBG%FHNDL,*) 'z0_t       =', z0_t
wwm_windinput.F90:158:                !WRITE(DBG%FHNDL,*) 'z_0        =', z0(ip)
wwm_windinput.F90:159:                !WRITE(DBG%FHNDL,*) 'Hs =', HS_w, 'L =', Lur, 'TP =', TP_w
wwm_windinput.F90:160:                !WRITE(DBG%FHNDL,*) 'Ulur       =' ,ULur, 'Ur =', Ur, 'EPS_D =', EPS_D
wwm_windinput.F90:161:                !WRITE(DBG%FHNDL,*) 'T_ds & T_t =', TAUW(ip), TAUHF(ip)
wwm_windinput.F90:162:                !WRITE(DBG%FHNDL,*) 'UFRIC      =', UFRIC1, UFRIC2
wwm_windinput.F90:246:              !WRITE(DBG%FHNDL,'(2I10,4F15.8)') IS, ID, SSINE(IS,ID), AUX3, AUX2, AUX1
