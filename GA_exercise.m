clear;
%遗传算法练习：最大值搜索
p = 250; %种群大小,必须为偶数
l = 30;%染色体尺寸(警告：该参数涉及进制转换问题，不可随意修改）
h = 2;%染色体数量
epoch = 150; %迭代次数限制
cal = 0;	  %迭代计数器
pc = 0.55; %选择比例
pm0 = 0.26;%初代基因变异率
pm1 = 0.03;%后代基因变异率
delta = 0.000000;%变异率衰减系数
pt = 0.99;%基因交换率
%迭代分析数据存储矩阵
Amax = zeros(epoch,1);

%搜索区域（以原点为中心的正方形）和初值
size =  4;
u = 4/power(2,l);
x0 = -2.25;
y0 = -2.25; 
%将数值转化为向量，去除二进制数的最后两位
ini_x = bitget(int32(x0/u),32:-1:3);
ini_y = bitget(int32(y0/u),32:-1:3);


%遗传信息存储矩阵
CHROM = zeros(p,l,h); %基因型存储矩阵
PHENO = zeros(p,h);	%表现型存储矩阵
%适应度向量
ADAPT = zeros(p,1);
%选择结果向量
CHOOSED = zeros(round(p*pc));
%繁殖矩阵
BREED = zeros(p,l,2);
%准备进制转换矩阵
V2 = ones(l,1);
for i = l:-1:2
    V2(i-1,1) = 2*V2(i,1);
end
V2(1,1) = -1*V2(1,1);%Matlab以补码表示负数，第一位权重取反



for i = 1:p	%装载初值
	CHROM(i,:,1) = ini_x;
	CHROM(i,:,2) = ini_y;
end

%遗传算法开始
while(cal<epoch)
    if cal == 0
        pm = pm0;
    else
        pm = pm1;
    end
    %更新基因变异率
  	pm = pm*exp(-1*delta*epoch);
    %计算染色体对应的表现型
    cal = cal+1;
    for i = 1:h;
    	PHENO(:,i) = CHROM(:,:,i)*V2; 
    end
    PHENO = PHENO*u*4;
    %计算适应度
	for i = 1:p
		ADAPT(i,1) = adaptability(PHENO(i,1),PHENO(i,2));
	end
	Amax(cal) = max(ADAPT);
	%轮盘选择
	if min(ADAPT)<0	%处理负适应度
    	ADAPT = ADAPT+abs(min(ADAPT));
	end
	s = sum(ADAPT);
	ADAPT = ADAPT/s;
	ADAPT = 1000*ADAPT;%建造轮盘
	for i = 2:p
    	ADAPT(i,1) = ADAPT(i-1,1)+ADAPT(i,1);
	end
	for i = 1:round(p*pc)%抽取次数
		Rvalue = rand()*1000;%获取轮盘值
		for j = 1:p	%查找对应的个体
    		if ADAPT(j) >= Rvalue
        		CHOOSED(i) = j;
  				break
     		end
		end
	end
	%复制被抽取的个体到繁殖矩阵
	for i = 1:round(p*pc)
    	BREED(i,:,:) = CHROM(CHOOSED(i),:,:);
	end
	for i = round(p*pc):p		%随机复制已被抽取的个体填补剩余的位置
    	j = ceil(rand()*p*pc);
    	BREED(i,:,:) = BREED(CHOOSED(j),:,:);
	end
	%基因交换
	for i = 1:(p/2)
		if rand()<=pt	%根据概率判断是否交换
 			for j = 1:h	%各个染色体单独交换
  	  		%Clength = fix(rand*l);		%抽取交换位点数量
  			Clength =  10;					%固定交换位点
			Cgen = BREED(2*i,(l-Clength):l,j) ;
 			BREED(2*i,(l-Clength):l,j) = BREED(2*i-1,(l-Clength):l,j);
 			BREED(2*i-1,(l-Clength):l,j) = Cgen;
 			clear Cgen; %删除保存交换值的临时变量
 			end
		end
	end
	%基因变异
	for i = 1:p			%个体
    	for j = 1:h		%染色体
        	for k = 1:l	%基因位点
             if rand()<=pm
                BREED(i,k,j) = ~BREED(i,k,j);
             end
        	end
    	end
	end
	CHROM = BREED;
end


%处理输出结果
fprintf('Analysing...\n');
plot(linspace(1,epoch,epoch),Amax);
xlabel('Epoch');
ylabel('Maximum Adaptability');