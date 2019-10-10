function [ x,result_CC,CC_time ] = function_CC( model,fid )
% Project Name(��Ŀ����
% ================
% Center-Cut algorithm
% ����Ŀ���ܣ�
% WHY?
% WHAT?
% �Ƽ�̼�ŷŽ��ף�carbon emission trading��CET��������ϣ�unit commitment��UC������ĸ�Ч��ⷽ���ǵ���ϵͳ�Ż������о��ȵ�֮һ��
% ���ĸ�ƽ�棨Center-Cut��CC)�㷨��Ϊһ�����Խ����㷨�����������Խ����㷨һ���蹹��������Լ�������Խ��ƶ����塣
% ����֮ͬ������CC�㷨Ѱ�Ҹ����Խ��ƶ������Chebyshev������Ϊ����⣬���ڸ������λ�ڶ������ڲ�����˼��п����ǿ��н⡣������һ����ʹ��CC�㷨�ܿ����ҵ����н⡣
% ���Ľ�CC�㷨�����������CPLEX���жԱȣ�ʵ����������CC�㷨��CPLEX�����ҵ����������н��Ҹ��ʺ������ģ����ϵͳ���⡣
% 
% In power system, it is one of the research hotspots that finding an efficient method for solving unit commitment (UC) problem considering carbon emission trading (CET). 
% As a linear approximation algorithm, Center-Cut (CC) algorithm also needs to construct approximate polyhedron of non-linear constraints like other linear approximation algorithms. 
% But the difference is that the algorithm chooses Chebyshev Center of the polyhedron as the trial solution. It is very likely to be a feasible solution, because the trial solution is interior point of the polyhedron. 
% In this paper, CC algorithm is compared with CPLEX, which is a commercial solver. 
% The experimental results show that CC algorithm can find high quality feasible solutions faster than CPLEX and is suitable for solving large-scale UC-CET problem.
% 
% User Guide���û�ָ�ϣ�
% -----------
% 
% Prerequisite(���и���Ŀ��Ҫ�Ŀ�������):
% -----------
% CPLEX
% matlab
% 
% Publication:
% -----------
% If you use Our ��Ŀ��/���� in academic work then please consider citing our papers.
% (������ο����ǵ���Ŀ/���� ������ѧ��Ŀ�ģ��뿼���������ǵ�����
% �������ĸ�ƽ���㷨���Ƽ�̼�ŷŽ��׵Ļ����������
% A Center-Cut Algorithm for Solving the Unit Commitment Problem Considering Carbon Emission Trading
% About Us
% -----------
% Authors����ά(Li Wei)  ���ַ�(Yang Linfeng)
% Team��
% Webpage: http://jians.gxu.edu.cn/default.do
TimeLimit=model.TimeLimit;
tolerances.TimeLimit=TimeLimit;
tolerances.solution_limit=0;
termination_LP=1e-3;
total_time=0;
iter_num=1;
termination_r=1e-3;
lambda=1;
L=model.L;

orgin_H=model.H;
orgin_f=model.f;
col_H=size(orgin_H,1);
[ model ] = obj_linear( model,L );
H=model.H;
f=[model.f;0];
Aineq=[model.Aineq,sparse(size(model.Aineq,1),1)];
bineq=model.bineq;
Aeq=[model.Aeq,sparse(size(model.Aeq,1),1)];
beq=model.beq;
l=[model.l;0];
Q=[model.Q,sparse(size(model.Q,1),1);sparse(1,size([model.Q,sparse(size(model.Q,1),1)],2))];
r=model.r;
lb=[model.lb;0];
ub=[model.ub;inf];
ctype=strcat(model.ctype,'C');
Value_num=length(ctype);
two_bin_location=strfind(ctype,'B');
NLP_Aineq=Aineq;
NLP_bineq=bineq;

ctype_NLP=ctype;
ctype_NLP(1,two_bin_location)='C';

result=0;
while(TimeLimit>0||iter_num<=2)%-1*result>termination_r
    if(iter_num~=1)
        if(gap_g<=0.001)
            Aeq_UC_init_test=sparse(1:length(two_bin_location),two_bin_location,1,length(two_bin_location ),Value_num);
            beq_UC_init_test=x(two_bin_location,1);
            tolerances.TimeLimit=TimeLimit;
            tolerances.TOLER=1e-3;
            [result,x,time ] = sovle( H,f,NLP_Aineq,NLP_bineq,[Aeq;Aeq_UC_init_test],[beq;beq_UC_init_test],Q,l,r,lb,ub,ctype_NLP,4,tolerances);
            if(isempty(x))
              x=x_k_1;
            total_time=total_time+time;
            CC_time=total_time;
            result_CC=obj;
                return;
            end
            deta_x=((x-x_k_1)'*(x-x_k_1))^0.5;
            x_k_1=x;
            total_time=total_time+time;
            TimeLimit=TimeLimit-time;
            gap_g=x'*Q*x+l'*x-r;
            obj=0.5*x(1:col_H,1)'*orgin_H*x(1:col_H,1)+orgin_f'*x(1:col_H,1);
            disp_result( fid,'CC_NLP',iter_num,gap_g,total_time,result,obj,deta_x);
            if(TimeLimit<0)
                CC_time=total_time;
                result_CC=obj;
                return;
            end
            %             gradient=sparse(x'*H+f');
            %             bineq_obj_cutting=gradient*x;
            %             gradient(1,Value_num)=lambda*(gradient*gradient')^0.5;
            %             Aineq_obj_cutting=gradient;
            if(isempty(H))
                gradient=f';
            else
                gradient=sparse(x'*H+f');
            end
            bineq_obj_cutting=gradient*x;
            gradient(1,Value_num)=lambda*(gradient*gradient')^0.5;
            Aineq_obj_cutting=gradient;
            
            Aineq_cutting=[];
            bineq_cutting=[];
            if(gap_g>-0.001&&gap_g<0.001)
                gradient=sparse(2*x'*Q+l');
                bineq_cutting=gradient*x-(x'*Q*x+l'*x-r);
                Aineq_cutting=gradient;
            end
            Aineq=[Aineq;Aineq_obj_cutting;Aineq_cutting];
            bineq=[bineq;bineq_obj_cutting;bineq_cutting];
        else
            gradient=sparse(2*x'*Q+l');
            bineq_cutting=gradient*x-(x'*Q*x+l'*x-r);
             gradient(1,Value_num)=(gradient*gradient')^0.5;
            Aineq_cutting=gradient;
            
            Aineq=[Aineq;Aineq_cutting];
            bineq=[bineq;bineq_cutting];
        end
        tolerances.TimeLimit=TimeLimit;
        tolerances.TOLER=5e-3;
        tolerances.solution_limit=2;
        [result,x,time ] = sovle( [],f_MILP,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,2,tolerances);
        if(isempty(x))
            x=x_k_1;
             total_time=total_time+time;
            CC_time=total_time;
            result_CC=obj;
            return;
        end
        deta_x=((x-x_k_1)'*(x-x_k_1))^0.5;
        x_k_1=x;
        
    else

        tolerances.TimeLimit=TimeLimit;
        tolerances.TOLER=5e-3;
        tolerances.solution_limit=2;
        [result,x,time ] = sovle( H,f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,2,tolerances);
        if(isempty(x))
            x=x_k_1;
            total_time=total_time+time;
            CC_time=total_time;
            result_CC=obj;
            return;
        end
        x_k_1=x;
        deta_x=((x-x_k_1)'*(x-x_k_1))^0.5;
        f_MILP=sparse(Value_num,1);
        f_MILP(Value_num,1)=-1;
    end
    gap_g=x'*Q*x+l'*x-r;
    obj=0.5*x(1:col_H,1)'*orgin_H*x(1:col_H,1)+orgin_f'*x(1:col_H,1);
    total_time=total_time+time;
    TimeLimit=TimeLimit-time;
    disp_result( fid,'CC_MILP',iter_num,gap_g,total_time,result,obj,deta_x);
    if(TimeLimit<0)
        CC_time=total_time;
        result_CC=obj;
        return;
    end
    iter_num=iter_num+1;
end
Aeq_UC_init_test=sparse(1:length(two_bin_location),two_bin_location,1,length(two_bin_location ),Value_num);
beq_UC_init_test=x(two_bin_location,1);
tolerances.TimeLimit=TimeLimit;
tolerances.TOLER=1e-3;
[result,x,time ] = sovle( H,f,NLP_Aineq,NLP_bineq,[Aeq;Aeq_UC_init_test],[beq;beq_UC_init_test],Q,l,r,lb,ub,ctype_NLP,4,tolerances);
total_time=total_time+time;
 deta_x=((x-x_k_1)'*(x-x_k_1))^0.5;
  gap_g=x'*Q*x+l'*x-r;
obj=0.5*x(1:col_H,1)'*orgin_H*x(1:col_H,1)+orgin_f'*x(1:col_H,1);
disp_result( fid,'CC_NLP',iter_num,gap_g,total_time,result,obj,deta_x);
result_CC=obj;
CC_time=total_time;

end

