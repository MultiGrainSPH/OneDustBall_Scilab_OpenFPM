pi=%pi


//данные задачи

T=0.2;
Rg0=1; //радиус газа
Mg=1000; //масса газа
Rd0=1; //радиус пыли
e0=18.6 // энергия (через нее можно вычислить ускорение газа)
//aRg0=10; //ускорение газа

VRg0=-1; //скорость газа
VRd0=-1; //скорость пыли не является свободным параметром при малых t_stop
Md=0; //масса пыли
ghamma=4.0/3.0;
t_stop=0.00001;





//разбиение
tau=0.000001; //T/tau должно быть целым, tau<t_stop
//if tau>=t_stop then
//    tau=t_stop*0.5
//    printf('Cлишком большое tau. Далее расчеты будут приведены при tau=%f', tau);
//end
//if int(T/tau)<T/tau then
//    tau=T/(int(T/tau)+1);
//    printf('T/tau не является целым. Далее расчеты будут приведены при tau=%f', tau);
//end
N=T/tau;

/////////////////////////////

//вычислим
aRd0=Rd0*(VRg0/Rg0-VRd0/Rd0)/t_stop;
otn_ro0=Md*Rg0**3/(Mg*Rd0**3)
aRg0 = 2*(ghamma-1)*e0/Rg0 - otn_ro0*Rg0*(VRg0/Rg0-VRd0/Rd0)/t_stop;
C=Rg0**(3*ghamma-1)*(aRg0/Rg0+otn_ro0*aRd0/Rd0)/(2*ghamma-2);

function f=f(ui)
    f=zeros(4,1);
    f(1) = ui(2);
    f(4) = ui(3)*(ui(2)/ui(1)-ui(4)/ui(3))/t_stop;
    ro_g=3*Mg/(4*pi*ui(1)**3);
    ro_d=3*Md/(4*pi*ui(3)**3);
    f(2) = 2*(ghamma-1)*C*ui(1)**(-3*ghamma+2)-ro_d*f(4)*ui(1)/(ro_g*ui(3));
    f(3) = ui(4);
 
endfunction

function u=Method_Euler(N)
    u=zeros(N+1, 4);
    u(1, 1)=Rg0;
    u(1, 2)=VRg0;
    u(1, 3)=Rd0;
    u(1, 4)=VRd0;
    //ti=0;
    for i=1:N;
        g=f(u(i,:))
        u(i+1, 1)=u(i, 1)+g(1)*tau;
        u(i+1, 2)=u(i, 2)+g(2)*tau;
        u(i+1, 3)=u(i, 3)+g(3)*tau;
        u(i+1, 4)=u(i, 4)+g(4)*tau;
        //ti=ti+tau;
    end
endfunction

u=Method_Euler(N);
t=[0:tau:T];
plot(t,u(:,1),'r');
plot(t,u(:,3));
//C0=aRg0*Rg0**(3*ghamma-2)/2;
//B = VRg0**2 + C0/(Rg0*Rg0);
//C1 = ((B*Rg0**2-C0)**0.5)/B;
//plot(t,(((t+C1)**2 * B**2 + C0)/B)**0.5, 'black');
xtitle   ("Красным Rg, синим Rd",   "t",   "R");

//вывод в файлы
M1=zeros(101, 7);
ii=0
for i=1:101
    M1(i,1)=t(ii+1);
    
    M1(i,2)=u(ii+1,1);
    M1(i,3)=u(ii+1,3);
    M1(i,4)=u(ii+1,2);
    M1(i,5)=u(ii+1,4);
    ii=ii+10000;
    M1(i,6)=3*Mg/(4*pi*u(i,1)**3);
    M1(i,7)=3*Md/(4*pi*u(i,3)**3);
end
//M1(:,1)=t;
//M1(:,2)=u(:,1);
//M1(:,3)=u(:,3);
//M1(:,4)=u(:,2);
//M1(:,5)=u(:,4);
//for i=1:N+1;
//    M1(:,6)=3*Mg/(4*pi*u(i,1)**3);
//    M1(:,7)=3*Md/(4*pi*u(i,3)**3);
//end
some_comment1 = ["t         Rg        Rd        Vg        Vd        ro_g      ro_d"];
fprintfMat("C:\Users\opskl\Downloads\File1.dat", M1, "%5.7f",some_comment1);
R=int(u(N+1,1));
if R<u(N+1,1) then
    R=R+1;
end
r=[0:R*0.01:R];
M2=zeros(101, 7)
for i=1:101;
    M2(i,1)=r(i);
    if r(i)<u(N+1,1) then
        M2(i,2)=3*Mg/(4*pi*u(N+1,1)**3);
        M2(i,3)=(ghamma-1)*M2(i,2)*(1-(r(i)/u(N+1,1))**2)*C*u(N+1,1)**(3-3*ghamma);
        M2(i,4)=(1-(r(i)/u(N+1,1))**2)*C*u(N+1,1)**(3-3*ghamma);
        M2(i,7)=u(N+1,2)*r(i)/u(N+1,1);
    else
        M2(i,2:4)=zeros(1,3);
        M2(i,7)=0;
    end
    if r(i)<u(N+1,3) then
        M2(i,5)=3*Md/(4*pi*u(N+1,3)**3);
        M2(i,6)=u(N+1,4)*r(i)/u(N+1,3);
    else
        M2(i,5:6)=zeros(1,2);
    end
end
some_comment2 = ["r         ro_g      p          e         ro_d       Vd        Vg"];
fprintfMat("C:\Users\opskl\Downloads\File2.dat", M2, "%5.7f",some_comment2);



