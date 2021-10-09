x=load('C:\Users\avboz\PycharmProjects\task_3b\x.txt');
y=load('C:\Users\avboz\PycharmProjects\task_3b\y.txt');
n=size(y, 2); % количество измерений в каждой точке

% #1.
% Вычислить в каждой точке xi
% 1) средние арифметические значения yi_
y_=zeros(length(x), 1);
for i = 1:length(x)
    y_(i, 1)=sum(y(i,:))./n;
end
writematrix(y_,'results.xlsx','Sheet','y_','Range','A1')
% 2) несмещенные оценки дисперсий si_2
s_2=zeros(length(x), 1);
for i = 1:length(x)
    s_2(i, 1)=(1/(n-1)).*sum((y(i,:)-y_(i,1)).^2);
end
writematrix(s_2,'results.xlsx','Sheet','s_2','Range','A1')
% 3) параметрические толерантные пределы для погрешностей
limits=zeros(length(x), 2);
% Толерантный множитель, взят из таблицы 4.2
% в учебнике "Таблицы математической статистики" Большев, Смирнов (стр. 238)
% gamma = Q = 0.95, P = 0.95, n = 10
k=3.379;
s=sqrt(s_2);
for i = 1:length(x)
    limits(i,1)=y_(i,1)-k*s(i,1);
    limits(i,2)=y_(i,1)+k*s(i,1);
end
% 4) доверительные интервалы для математических ожиданий
interval=zeros(length(x), 2);
Q=0.95;
% Коэффициент Стъюдента
t=tinv((1+Q)/2, n-1);
for i = 1:length(x)
    interval(i,1)=y_(i,1)-(s(i,1)/sqrt(n))*t;
    interval(i,2)=y_(i,1)+(s(i,1)/sqrt(n))*t;
end
% Проверить гипотезу о равенстве дисперсий в этих точках по критерию Кочрена
s_2_max=0;
for i = 1:length(x)
    if s_2(i,1)>s_2_max
        s_2_max=s_2(i,1);
    end
end
g=s_2_max/sum(s_2);
% Критическое значение g_alpha(n-1=10-1=9, k=size(x)=41), alpha=0.95
% выбирается из таблиц критических значений критерия Кочрена
g_alpha=0.0745;
if g <= g_alpha
    str=sprintf('%.4f <= %.4f',g,g_alpha);
    disp(str);
    disp('Нулевая гипотеза НЕ противоречит экспериментальным данным, и для аппроксимации применяется МНК')
else
    str=sprintf('%.4f > %.4f',g,g_alpha);
    disp(str);
    disp('Нулевая гипотеза противоречит экспериментальным данным, и для аппроксимации применяется ОМНК')
end

% #2.
% Произвести последовательную полиномиальную аппроксимацию.
% В качестве значений y при аппроксимации необходимо использовать средние арифметические значения y_.
k = length(x); % число точек, в которых проводились измерения
alpha=0.05;
flag = true;
% #2.1. 
% Начать с нулевой степени полинома
% #2.2. 
% Вычислить оценки коэффициентов полинома a 
% с помощью МНК для заданной степени полинома
a=zeros(k, k);
condition=zeros(k, 1);
F=zeros(k-1, 2);
F_check=strings(k-1, 1);
for q = 0:k-2
    str=sprintf('q = %d',q);
    disp(str);
    % Формируется матрица X
    X=zeros(k, q+1);
    for i = 0:q
        X(:,i+1)=x.^i;
    end
    % Оцениваются коэффициенты полинома a
    X_T=transpose(X);
    mul=X_T*X;
    % coeff=0.0;
    % E=diag(diag(ones(q+1)));
    % add=mul+coeff*E;
    condition(q+1,1)=cond(mul);
    inverted=inv(mul);
    mul2=inverted*X_T;
    a(1:q+1,q+1)=mul2*y_;
    % #2.3.
    % Проверить гипотезу о степени q полинома
    s_epsilon_2=sum(s_2)/k;
    mul3=X*a(1:q+1,q+1);
    sub=mul3-y_;
    sub_T=transpose(sub);
    mul4=sub_T*sub;
    R_2=(n/s_epsilon_2)*mul4;
    F(q+1,1)=R_2/(k-q-1);
    % Для проверки гипотезы используется критерий Фишера
    F(q+1,2)=finv(1-alpha, k-q-1, n-1);
    if F(q+1,1) > F(q+1,2)
        F_check(q+1,1)='Неверна';
    else
        F_check(q+1,1)='Верна';
        % если гипотезу о степени q полинома не будет отвергнута
        % вычисляется ковариационная матрица
        if flag
            flag=false;
            % оценить ковариационную матрицу оценок коэффициентов S_a и дисперсии s_ai_2 
            % (это значения на главной диагонали ковариационной матрицы)
            S_a_q=(s_epsilon_2/n)*inverted;
            writematrix(S_a_q,'results.xlsx','Sheet','S_a_q','Range','A1');
            % #2.4.
            % Вычислить корреляционную матрицу R_a и
            % коэффициенты корреляции r_a_ij между оценками коэффициентов по матрице ковариации
            % Корреляционная матрица — это матрица, элементы которой являются коэффициентами корреляции
            R_a_q=zeros(q+1,q+1);
            for i = 1:q+1
                for j = 1:q+1
                    R_a_q(i,j)=S_a_q(i,j)/sqrt(S_a_q(i,i)*S_a_q(j,j));
                end
            end
            writematrix(R_a_q,'results.xlsx','Sheet','R_a_q','Range','A1');
            condition_q=cond(mul);
            % #2.5.
            % Была получена степень q полинома, прошедшая гипотезу о степени полинома. 
            % Произвести все те же действия для полинома степени, равной k-1 
            % (вычислить коэффициенты и корреляцию между ними)
            str=sprintf('q = k-1 = %d',k-1);
            disp(str);
            % Формируется матрица X
            X=zeros(k, k);
            for i = 0:k-1
                X(:,i+1)=x.^i;
            end
            % Оцениваются коэффициенты полинома a
            X_T=transpose(X);
            mul=X_T*X;
            coeff=0.0;
            E=diag(diag(ones(k)));
            add=mul+coeff*E;
            inverted=inv(add);
            mul2=inverted*X_T;
            a(1:k,k)=mul2*y_;
            S_a_k_1=(s_epsilon_2/n)*inverted;
            writematrix(S_a_k_1,'results.xlsx','Sheet','S_a_k_1','Range','A1');
            R_a_k_1=zeros(k,k);
            for i = 1:k
                for j = 1:k
                    R_a_k_1(i,j)=S_a_k_1(i,j)/sqrt(S_a_k_1(i,i)*S_a_k_1(j,j));
                end
            end
            writematrix(R_a_k_1,'results.xlsx','Sheet','R_a_k_1','Range','A1');
            condition(k,1)=cond(add);
            condition_k_1=cond(add);
        end
        figure;
        p=0;
        for i = 0:q
            p=p+a(i+1,q+1)*x.^i;
        end
        hold on;
        p1 = plot(x, p, 'Color','m');
        p2 = plot(x, y_, 'o', 'MarkerEdgeColor','r', 'MarkerSize',4);
        for i = 1:10
            p3 = plot(x, y(:,i), 'o', 'MarkerEdgeColor','g', 'MarkerSize',2);
        end
        p4 = plot(x, limits(:,1), 'Color','blue', 'LineWidth',1.1);
        plot(x, limits(:,2), 'Color','blue', 'LineWidth',1.1);
        p5 = plot(x, interval(:,1), 'Color','black', 'LineWidth',1.1);
        plot(x, interval(:,2), 'Color','black', 'LineWidth',1.1);
        p=0;
        for i = 0:k-1
            p=p+a(i+1,k)*x.^i;
        end
        p6 = plot(x, p, 'Color','c');
        s1 = sprintf('Полином %d-ой (q-ой) степени', q);
        s2 = sprintf('Полином %d-ой ((k-1)-ой) степени', k-1);
        legend([p1 p2 p3 p4 p5 p6],{s1,'Точки, по которым проводилась аппроксимация','Исходные точки','Толерантные пределы','Доверительные интервалы',s2});
        ylim([-250 450]);
        ax = gca;
        str = sprintf('2.%d.jpg', q);
        exportgraphics(ax,str);
        close;
    end
end
writematrix(condition,'results.xlsx','Sheet','cond','Range','A1');
writematrix(a,'results.xlsx','Sheet','a','Range','B2');
writematrix(F,'results.xlsx','Sheet','F','Range','A1');
writematrix(F_check,'results.xlsx','Sheet','F','Range','C1');

% #3.
% Произвести аппроксимацию исходной зависимости другими способами
% Представить полученные графики аппроксимации
% #3.1.
% Произвести аппроксимацию зависимости прямой линией с помощью функций
% 1) regress (использует метод R-Square)
[b,~]=regress(y_,X(:,1:2),0.05);
Y=b(2)*x+b(1);
figure;
hold on;
p1 = plot(x, Y, 'Color','c', 'LineWidth',1.1);
% 2) robustfit (робастная регрессия)
b=robustfit(x,y_);
Y=b(2)*x+b(1);
p2 = plot(x, Y, 'Color','b', 'LineWidth',1.1);
% 3) polyfit (полиномиальная регрессия с n=1)
n=1;
b=polyfit(x,y_,n);
Y=b(1)*x+b(2);
p3 = plot(x, Y, '--', 'Color','#A2142F', 'LineWidth',1.1);
% 4) ridge (ридж-регрессия с регуляризацией)
b=ridge(y_,X(:,1:2),0.05);
Y=b(2)*x+b(1);
p4 = plot(x, Y, '--', 'Color','#EDB120', 'LineWidth',1.1);
% Представить полученные графики аппроксимации 
% (полученная аппроксимирующая кривая одним цветом, точки, 
% по которым проводилась аппроксимация маркерами одного типа 
% и все исходные точки маркерами другого типа)
p5 = plot(x, y_, 'o', 'MarkerEdgeColor','r', 'MarkerSize',4);
for i = 1:10
    p6 = plot(x, y(:,i), 'o', 'MarkerEdgeColor','g', 'MarkerSize',2);
end
legend([p1 p2 p3 p4 p5 p6],{'regress','robustfit','polyfit (n=1)','ridge','точки, по которым проводилась аппроксимация','исходные точки'});
ylim([-200 300]);
ax = gca;
exportgraphics(ax,'3.1.jpg');
close;

% #3.2.
% Произвести полиномиальную аппроксимацию с помощью функций polyfit (polyval). 
% Можно воспользоваться утилитой polytool, являющейся графическим интерфейсом к polyfit. 
% Подобрать степень полинома, наилучшим способом аппроксимирующую исходную зависимость.
n=19;
b=polyfit(x,y_,n);
figure;
hold on;
p1 = plot(x, polyval(b, x), 'Color','#A2142F', 'LineWidth',1.1);
p2 = plot(x, y_, 'o', 'MarkerEdgeColor','r', 'MarkerSize',4);
for i = 1:10
    p3 = plot(x, y(:,i), 'o', 'MarkerEdgeColor','g', 'MarkerSize',2);
end
legend([p1 p2 p3],{'polyfit (n='+string(n)+')','точки, по которым проводилась аппроксимация','исходные точки'});
ylim([-200 300]);
ax = gca;
exportgraphics(ax,'3.2.jpg');
close;
alpha=0.05;
n=19;
polytool(x,y_,n,alpha);

% #3.3.
% Произвести кусочную полиномиальную аппроксимацию с помощью функций
figure;
xq = -2.0:0.000001:2.0;
% 1) interp1 (линейная)
vq1=interp1(x,y_,xq,'linear');
hold on;
p1 = plot(xq,vq1, 'Color','c', 'LineWidth',1.1);
% 2) interp1 (кубическая)
vq2=interp1(x,y_,xq,'cubic');
p2 = plot(xq,vq2, 'Color','b', 'LineWidth',1.1);
% 3) pchip (полиномами Эрмита)
vq3=interp1(x,y_,xq,'pchip');
p3 = plot(xq,vq3, 'Color','#A2142F', 'LineWidth',1.1);
% 4) spline (сплайны)
vq4=interp1(x,y_,xq,'spline');
p4 = plot(xq,vq4, 'Color','#EDB120', 'LineWidth',1.1);
p5 = plot(x, y_, 'o', 'MarkerEdgeColor','r', 'MarkerSize',4);
for i = 1:10
    p6 = plot(x, y(:,i), 'o', 'MarkerEdgeColor','g', 'MarkerSize',2);
end
legend([p1 p2 p3 p4 p5 p6],{'linear','cubic','pchip','spline','точки, по которым проводилась аппроксимация','исходные точки'});
ylim([-200 300]);
xlim([-2.0 2.0]);
ax = gca;
exportgraphics(ax,'3.3.jpg');
close;
% #3.4.
% Произвести нелинейную аппроксимацию с помощью функции nlinfit. 
% В качестве нелинейной функции использовать 
% произведение полинома на гармоническую функцию:
% y = (sin(alpha*x+beta))*(a_n*x.^n+...+a_1*x+a_0)
alpha_beta=zeros(23, 3);
figure;
hold on;
modelfun = @(coeffs, x) (sin(coeffs(1)*x+coeffs(2)).*(coeffs(3)));
initials = ones(3,1);
new_coeffs = nlinfit(x, reshape(y_,[1,length(y_)]), modelfun, initials);
alpha_beta(1:3,1)=new_coeffs;
Y=modelfun(new_coeffs, x);
p1 = plot(x, Y, 'Color',[0.83 0.14 0.14], 'LineWidth',1.1);
modelfun = @(coeffs, x) (sin(coeffs(1)*x+coeffs(2)).*(coeffs(3)*x.^10+coeffs(4)*x.^9+coeffs(5)*x.^8+coeffs(6)*x.^7+coeffs(7)*x.^6+coeffs(8)*x.^5+coeffs(9)*x.^4+coeffs(10)*x.^3+coeffs(11)*x.^2+coeffs(12)*x+coeffs(13)));
initials = ones(13,1);
new_coeffs = nlinfit(x, reshape(y_,[1,length(y_)]), modelfun, initials);
alpha_beta(1:13,2)=new_coeffs;
Y=modelfun(new_coeffs, x);
p2 = plot(x, Y, 'Color',[1.00 0.54 0.00], 'LineWidth',1.1);
modelfun = @(coeffs, x) (sin(coeffs(1)*x+coeffs(2)).*(coeffs(3)*x.^20+coeffs(4)*x.^19+coeffs(5)*x.^18+coeffs(6)*x.^17+coeffs(7)*x.^16+coeffs(8)*x.^15+coeffs(9)*x.^14+coeffs(10)*x.^13+coeffs(11)*x.^12+coeffs(12)*x.^11+coeffs(13)*x.^10+coeffs(14)*x.^9+coeffs(15)*x.^8+coeffs(16)*x.^7+coeffs(17)*x.^6+coeffs(18)*x.^5+coeffs(19)*x.^4+coeffs(20)*x.^3+coeffs(21)*x.^2+coeffs(22)*x+coeffs(23)));
initials = ones(23,1);
new_coeffs = nlinfit(x, reshape(y_,[1,length(y_)]), modelfun, initials);
alpha_beta(1:23,3)=new_coeffs;
Y=modelfun(new_coeffs, x);
p3 = plot(x, Y, 'Color','c', 'LineWidth',1.1);
p4 = plot(x, y_, 'o', 'MarkerEdgeColor','r', 'MarkerSize',4);
for i = 1:10
    p5 = plot(x, y(:,i), 'o', 'MarkerEdgeColor','g', 'MarkerSize',2);
end
legend([p1, p2, p3, p4, p5],{'nlinfit n=0','nlinfit n=10','nlinfit n=20','точки, по которым проводилась аппроксимация','исходные точки'});
ylim([-200 300]);
xlim([-2.0 2.0]);
ax = gca;
exportgraphics(ax,'3.4.jpg');
close;