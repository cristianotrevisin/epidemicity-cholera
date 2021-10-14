close all
clear all
clc

load data/geodata POPnodes WS_dept
POPnodes = POPnodes' * WS_dept;
POP = POPnodes;

[rSeq1, rSeq2] = load_data();
opt = "best";
[x1,x2] = get50(opt,rSeq1, rSeq2);

[cases_AD1_week, time, y] = SIARBV(2,2,x1,x2);


figure()
subplot(3,4,1)
hold on
for i = 1:10
    plot(y(:,1+14*(i-1)));
end
title('S0')

subplot(3,4,2)
hold on
for i = 1:10
    plot(y(:,2+14*(i-1)));
end
title('I0')

subplot(3,4,3)
hold on
for i = 1:10
    plot(y(:,3+14*(i-1)));
end
title('A0')

subplot(3,4,4)
hold on
for i = 1:10
    plot(y(:,4+14*(i-1)));
end
title('R0')

subplot(3,4,5)
hold on
for i = 1:10
    plot(y(:,7+14*(i-1)));
end
title('S1')

subplot(3,4,6)
hold on
for i = 1:10
    plot(y(:,8+14*(i-1)));
end
title('I1')

subplot(3,4,7)
hold on
for i = 1:10
    plot(y(:,9+14*(i-1)));
end
title('A1')

subplot(3,4,8)
hold on
for i = 1:10
    plot(y(:,10+14*(i-1)));
end
title('R1')

subplot(3,4,9)
hold on
for i = 1:10
    plot(y(:,11+14*(i-1)));
end
title('S2')

subplot(3,4,10)
hold on
for i = 1:10
    plot(y(:,12+14*(i-1)));
end
title('I2')

subplot(3,4,11)
hold on
for i = 1:10
    plot(y(:,13+14*(i-1)));
end
title('A2')

subplot(3,4,12)
hold on
for i = 1:10
    plot(y(:,14+14*(i-1)));
end
title('R2')




load data/geodata POPnodes WS_dept
POPnodes = POPnodes' * WS_dept;
POP = POPnodes';


figure()
subplot(3,4,1)
hold on
for i = 1:10
    plot(y(:,1+14*(i-1))/POP(i));
end
title('S0')

subplot(3,4,2)
hold on
for i = 1:10
    plot(y(:,2+14*(i-1))/POP(i));
end
title('I0')

subplot(3,4,3)
hold on
for i = 1:10
    plot(y(:,3+14*(i-1))/POP(i));
end
title('A0')

subplot(3,4,4)
hold on
for i = 1:10
    plot(y(:,4+14*(i-1))/POP(i));
end
title('R0')

subplot(3,4,5)
hold on
for i = 1:10
    plot(y(:,7+14*(i-1))/POP(i));
end
title('S1')

subplot(3,4,6)
hold on
for i = 1:10
    plot(y(:,8+14*(i-1))/POP(i));
end
title('I1')

subplot(3,4,7)
hold on
for i = 1:10
    plot(y(:,9+14*(i-1))/POP(i));
end
title('A1')

subplot(3,4,8)
hold on
for i = 1:10
    plot(y(:,10+14*(i-1))/POP(i));
end
title('R1')

subplot(3,4,9)
hold on
for i = 1:10
    plot(y(:,11+14*(i-1))/POP(i));
end
title('S2')

subplot(3,4,10)
hold on
for i = 1:10
    plot(y(:,12+14*(i-1))/POP(i));
end
title('I2')

subplot(3,4,11)
hold on
for i = 1:10
    plot(y(:,13+14*(i-1))/POP(i));
end
title('A2')

subplot(3,4,12)
hold on
for i = 1:10
    plot(y(:,14+14*(i-1))/POP(i));
end
title('R2')

figure()
hold on
for i = 1:10
    plot(y(:,5+14*(i-1)));
end
title('B')