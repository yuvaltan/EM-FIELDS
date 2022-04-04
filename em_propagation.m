close all
clear all
clc
global delta h_out w_out w_in d_in H_out W_out W_in D_in h1 h2 
global iterations avg_runtime itration_avg_time_vector num_of_iteration_vector
delta = 0.1*(10^-3);
h_out = 4*(10^-3);
w_out = 10*(10^-3);
w_in = 0.8*(10^-3);
d_in = 2*(10^-3);
h1 = 0;
h2 = 0;
%% 
Set_values_by_delta(delta);

%--------------------------
% potential on the board
%--------------------------
board = board_line_potential();
board = all_board_potential(board);
number_of_iterations = iterations
   figure 
   grid on
   imagesc(board)    %plot potential
   gamma1 = colorbar;
   gamma1.Label.String = 'Potential';
   x = linspace(0,w_out/delta,5); 
   y = linspace(0,h_out/delta,10);
   xticks(x)
   yticks(y)
   xlabel('y [0.1*mm]')
   ylabel('x [0.1*mm]')
   title('board potential')
   

%-------------------------- 
% electric field
%--------------------------
wires_area1 = wires_area();          
electric_field = the_electric_field(board,wires_area)
   figure 
   quiver(real(electric_field),imag(electric_field),5);  %plot electric field
   title('electric field')       
   xticks(x)
   yticks(y)
   xlabel('y [0.1*mm]')
   ylabel('x [0.1*mm]')

%-------------------------- 
% magnetic field
%-------------------------- 
gamma = gamma_by_integral(electric_field);
magnetic_field1 = the_magnetic_field(electric_field,gamma,wires_area);
   figure
   quiver(real(magnetic_field1),imag(magnetic_field1),5) %plot magnetic field
   title('magnetic field')
   xticks(x)
   yticks(y)
   xlabel('y [0.1*mm]')
   ylabel('x [0.1*mm]')
   
%--------------------------------------------- 
% The influence of h1 h2 on gamma
%---------------------------------------------  
gamma_offset_mat=gamma_matrix();
   figure
   imagesc(gamma_offset_mat)
   gamma1 = colorbar;
   gamma1.Label.String = 'gamma';
   title('gamma(h1,h2)')
   x = linspace(0,5,5); 
   y = linspace(0,5,5);
   xticks(x)
   yticks(y)
   xlabel('h1')
   ylabel('h2')

%----------------------------------------------------------
% C, L, Zc, v (for regular board, and the changed board) 
%----------------------------------------------------------
h1=6;
h2=0;
gamma1=gamma;
Set_values_by_delta(0.1*(10^-3));
board = board_line_potential();
board = all_board_potential(board);
wires_area1 = wires_area();          
electric_field = the_electric_field(board,wires_area)
gamma6 = gamma_by_integral(electric_field);
u = 4*pi*1e-7;
epsilon = 2.5*1e-9/(36*pi) ;
eta = sqrt(u/epsilon);
%[C,L,Zc,v]
line_parameters1 = [epsilon*gamma1,u/gamma1,eta/gamma1,1/sqrt(u*epsilon)];
line_parameters6 = [epsilon*gamma6,u/gamma6,eta/gamma6,1/sqrt(u*epsilon)];

%--------------------------------------------- 
% hedim chart 
%--------------------------------------------- 
Zc1 = line_parameters1(3);
Zc6 = line_parameters6(3);
T=50;
hedim_mat = hedim_chart(Zc1,Zc6,T);
   figure
   imagesc(hedim_mat) 
   gamma1 = colorbar;
   gamma1.Label.String = '20Log|V/1[v]| [dB]';
   xlabel('z')
   ylabel('time')
   title('hedim chart [dB]')




%functions

%-----------------------------------------------
% defy sizes - for all modifications of delta
%-----------------------------------------------
function Set_values_by_delta(delta1)
global H_out W_out W_in D_in h_out w_out w_in d_in
H_out = h_out/delta1; 
W_out = w_out/delta1; 
W_in = w_in/delta1;
D_in = d_in/delta1;
end

%--------------------------
% lines potential
%--------------------------
function matrix = board_line_potential()
global  H_out W_out W_in D_in h1 h2
board = zeros (H_out+2*h2,W_out+2*h2);
for i = (W_out/2)-(D_in/2)-W_in+h2-h1:(W_out/2)-(D_in/2)+h2+h1   %a to b
    for j = (H_out/2)-(W_in/2)+h2-h1:(H_out/2)+(W_in/2)+h2+h1    %f to g
         board (j,i) =-1/2; 
    end
end
for i = (W_out/2)+(D_in/2)+h2-h1:(W_out/2)+(D_in/2)+W_in+h2+h1    %c to d
       for j = (H_out/2)-(W_in/2)+h2-h1:(H_out/2)+(W_in/2)+h2+h1  %f to g
          board (j,i) =1/2; 
       end
end
matrix = board;
end

%-----------------
% wires area
%-----------------
function wires_area_arr = wires_area()
    global H_out W_out W_in D_in h1 h2
    arr = zeros(2,(W_in+2*h1+1)^2);    %save the wire points
    k=0;
    for im = (H_out/2-W_in/2+h2-h1):(H_out/2+W_in/2+h2+h1)            %f to g
        for re = (W_out/2-D_in/2-W_in+h2-h1):(W_out/2-D_in/2+h2+h1)     %a to b
            k=k+1;
            arr(1,k) = re+im*(1i);                                
        end
    end
    k=0;
    for im = (H_out/2-W_in/2+h2-h1):(H_out/2+W_in/2+h2+h1)            %f to g
        for re=(W_out/2+D_in/2+h2-h1):(W_out/2+D_in/2+W_in+h2+h1)       %c to d
            k=k+1;
            arr(2,k) = re+im*(1i);                                
        end
    end
    wires_area_arr=arr;
end

%---------------------------------------------  
% the whole potential by numeric laplace
%--------------------------------------------- 
function potential = all_board_potential(board)
global  H_out W_out W_in D_in normVal h1 h2 iterations avg_runtime
normVal = inf;
mat1 = board; 
num_iter = 0; sum_iter_time=0;
while true
    tic;
    num_iter=num_iter+1;
    board = mat1;
    for j = 2:W_out+2*h2-1        %0 to wout
        for i = 2:(H_out/2-W_in/2+h2-h1)-1         %under the wires
            mat1(i,j) = (board(i-1,j)+board(i+1,j)+board(i,j-1)+board(i,j+1))/4;    %numeric laplace
        end
    end
    for j = 2:W_out+2*h2-1       
        for i = H_out/2+W_in/2+h2+h1+1:H_out+2*h2-1 %above the wires
            mat1(i,j) = (board(i-1,j)+board(i+1,j)+board(i,j-1)+board(i,j+1))/4;     
        end
    end
    for i = (H_out/2-W_in/2-1+h2-h1):(H_out/2+W_in/2+h2+h1+1)      %f to g
        for j = 2:(W_out/2-D_in/2-W_in)-1+h2-h1                    %start to a
            mat1(i,j) = (board(i-1,j)+board(i+1,j)+board(i,j-1)+board(i,j+1))/4;    
        end
        for k = (W_out/2-D_in/2+1+h2+h1):(W_out/2+D_in/2+h2-h1-1)     %b to c
            mat1(i,k) = (board(i-1,k)+board(i+1,k)+board(i,k-1)+board(i,k+1))/4;    
        end
        for m = (W_out/2+D_in/2+W_in+h2+h1+1):W_out+2*h2-1        %d to end
            mat1(i,m) = (board(i-1,m)+board(i+1,m)+board(i,m-1)+board(i,m+1))/4;    
        end
        for k = (W_out/2-D_in/2-W_in+h2-h1):(W_out/2-D_in/2+h2+h1)   % left wire
            if i >=H_out/2-W_in/2+h2-h1     %down left lip
                break
            end
            mat1(i,k) = (board(i-1,k)+board(i+1,k)+board(i,k-1)+board(i,k+1))/4;    
        end
        for k = (W_out/2+D_in/2+h2-h1):(W_out/2+D_in/2+W_in+h2+h1)    % right wire
            if i <= H_out/2+W_in/2+h2+h1     %up lip
                break
            end
            mat1(i,k) = (board(i-1,k)+board(i+1,k)+board(i,k-1)+board(i,k+1))/4;    
        end
    end
    t1=toc;
    normVal = norm(mat1-board)/norm(board); 
    sum_iter_time = sum_iter_time+t1;
    if normVal < 1e-5    %until error<10_-5
        break
    end
end
iterations = num_iter;    %counting iterations
avg_runtime = sum_iter_time/num_iter;    %avg runtime for one iteration
potential = mat1;
end

%-------------------------
% electric field
%-------------------------
function matrix = the_electric_field(board,wires_area)
global delta   H_out W_out h2 
electric_field_mat = zeros(H_out+2*h2,W_out+2*h2); 
    for j = 2:W_out+2*h2-1         
        for i = 2:H_out+2*h2-1     
            index_ij = i*1i+j;
            if ismember(index_ij,wires_area)
                continue
            end    
            electric_field_mat(i,j) = ((board(i,j-1)-board(i,j+1))/(2*delta)+1i*(board(i-1,j)-board(i+1,j))/(2*delta)) ;    %numeric div
        end
    end
    matrix =electric_field_mat;
end

%-------------
% Gamma
%-------------
function gamma1 = gamma_by_integral(electric_field)
global W_out  D_in W_in  H_out delta h1 h2
bottom = 0; high = 0; left = 0; right = 0;

for i = (W_out/2)-(D_in/2)-W_in-1+h2-h1:(W_out/2)-(D_in/2)+h2+h1 %a to b bottom
    bottom=bottom+1;
    bot_vec(bottom) = imag(electric_field ((H_out/2)-(W_in/2)+h2-h1-1,i));
end
for i = (W_out/2)-(D_in/2)-W_in+h2-h1:(W_out/2)-(D_in/2)+h2+h1+1 %a to b top 
    high=high+1;
    high_vec(high) =-imag(electric_field ((H_out/2)+(W_in/2)+h2+h1+1,i));
end

for j = (H_out/2)-(W_in/2)+h2-h1 :(H_out/2)+(W_in/2)+h2+h1 +1    %f to g left
    left=left+1;
    left_vec(left) = real(electric_field(j,(W_out/2)-(D_in/2)-W_in+h2-h1-1));
end
for j = (H_out/2)-(W_in/2)+h2-h1 -1:(H_out/2)+(W_in/2)+h2+h1     %f to g right
    right=right+1;
    right_vec(right) =-real(electric_field(j,(W_out/2)-(D_in/2)+h2+h1+1));
end        
    gamma1 =(delta)*((trapz(bot_vec)+ trapz(high_vec)+ trapz(left_vec)+ trapz(right_vec))); 
end

%---------------------------------------------  
% magnetic field
%--------------------------------------------- 
function matrix = the_magnetic_field(electric_field,gamma,wires_area)
global H_out W_out h2 
magnetic_field_mat = zeros(H_out+2*h2,W_out+2*h2); 
for j = 2:W_out+2*h2-1         
        for i = 2:H_out+2*h2-1     
            index_ij = i*1i+j;
            if ismember(index_ij,wires_area)
                continue
            end    
            magnetic_field_mat(i,j) = 1/gamma * (-imag(electric_field(i,j))+1i*real(electric_field(i,j)));   
        end
end
    
    matrix = magnetic_field_mat;
end

%---------------------------------- 
% different gammas matrix
%---------------------------------- 
function gamma_mat = gamma_matrix()
global h1 h2
h_vector = [0,1,2,4];
gamma = zeros(4,4);
for j=1:4
    h1=h_vector(j);
    for i=1:4
        h2=h_vector(i);
        board = board_line_potential();
        board = all_board_potential(board);
        wires_area1 = wires_area();
        modal_electric_field = the_electric_field(board,wires_area1); % for every h, calculate gamma
        gamma(j,i) = gamma_by_integral(modal_electric_field);
    end
end
gamma_mat = gamma;
end

%----------------------
% hedim chart maker
%----------------------
function matrix = hedim_chart(Zc1,Zc2,T)
pulse = T/10;
hedim = zeros(8*T,1000);  %z(t)
Rg = Zc1*(313361602/(313361602+312463938)); %ID formula
D=Zc1/(Zc1+Rg);   %voltage divider 
return12 = (Zc2-Zc1)/(Zc1+Zc2); % return factor12
transfer12 = 1+return12;     %transfer factor12
source_return = (Rg-Zc1)/(Rg+Zc1);  %return from left (source)
pulseLeft = 1-pulse; 
pulseRight = 1;
backPulseLeft = 500;
backPulseRight = 500+pulse;
Vg = D; 
for t=1:(8*T)
    if mod(t,100)==1 && t~=1  %every time he starts move right
        pulseLeft = 1-pulse;
        pulseRight = 1;
        backPulseLeft =500;
        backPulseRight =500+pulse; 
        Vg = Vg*source_return*return12;
    end
    for z=pulseLeft:pulseRight     %begin pulse to finish
        if z<=0 || z>1000
            continue
        end
        hedim(t,z) = Vg; 
        if mod(t,2*T)>=T+1 && z>500      %if he is on the way back
           hedim(t,z) = hedim(t,z)*transfer12; 
        end
    end
    pulseLeft = pulseLeft+(500/T); 
    pulseRight = pulseRight +(500/T);
    if mod(t,2*T)<=T    %if he is on the way forward
            continue
    end
    for z=backPulseLeft:backPulseRight   %begin back pulse to finish
        if z<=0 || z>500
            continue
        end
        hedim(t,z) = hedim(t,z) + Vg*return12;
    end
    backPulseLeft=backPulseLeft-(500/T);
    backPulseRight=backPulseRight-(500/T);
end
matrix = 20*log10(abs(hedim));  %log chart
end