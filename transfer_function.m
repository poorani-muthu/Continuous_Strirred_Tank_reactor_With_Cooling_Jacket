
systf=tf(reac_cont) %convert to a 2x2 transfer function

figure(1)
sgtitle("Step Response of invidual Outputs wrt individual Inputs")

subplot(2,2,1)
step(systf(1,1)) % compute and plot step response of output 1 to input 1
title("step response of output 1 to input 1")
ax = gca;
ax.FontSize = 4;
subplot(2,2,2)
step(systf(2,1)) % compute and plot step response of output 2 to input 1
title("step response of output 2 to input 1")
subplot(2,2,3)
step(systf(1,2))
title("step response of output 1 to input 2")
subplot(2,2,4)
step(systf(2,2))
title("step response of output 2 to input 2")


C0 = pidstd(1,1,1); 
C = pidtune(systf(2,2),C0)