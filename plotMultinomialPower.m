numSample = 10000;
alpha = 0.05;

subplot(2, 3, 1);
multinomialPower(50, 2, alpha, numSample);
grid on;
title('m = 50, n = 2');

subplot(2, 3, 2);
multinomialPower(50, 5, alpha, numSample);
grid on;
title('m = 50, n = 5');

subplot(2, 3, 3);
multinomialPower(50, 10, alpha, numSample);
grid on;
title('m = 50, n = 10');

subplot(2, 3, 4);
multinomialPower(500, 20, alpha, numSample);
grid on;
title('m = 250, n = 10');

subplot(2, 3, 5);
multinomialPower(500, 50, alpha, numSample);
grid on;
title('m = 250, n = 25');

subplot(2, 3, 6);
multinomialPower(500, 100, alpha, numSample);
grid on;
title('m = 250, n = 50');