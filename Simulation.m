clear
close all
clc
format long;

load('trueRandomNumbers');

c = clock;
% time variable can be used for random seed value
time = round(10000*(c(6)-floor(c(6))));
clear c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Middle Square Generator
% parameters: seed, hist, interval[0 1]
% returns generated random numbers

% MSrandomNumbers = middleSquare(761462387, 1, 0);
% bitmapGenerator(100, MSrandomNumbers, 0);
% frequencyMonobitTest(MSrandomNumbers);
% runsTest(MSrandomNumbers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear Congruential Generator
% parameters: m, a, c, seed, number of numbers to generate, histogram, interval[0, 1]
% returns generated random numbers

% LCrandomNumbers = linearCongruential(2^32, 214013, 2531011, time, 1e4, 1, 1);
% bitmapGenerator(100, LCrandomNumbers, 0);
% frequencyMonobitTest(LCrandomNumbers);
% runsTest(LCrandomNumbers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inverse Transform Sampling - Exponential
% parameters: number of numbers to generate, number of delay time, interval[0, 1]
% returns generated random numbers

% ITExpRandomNumbers = inverseTransformExponential(1e2, 0, 1);
% bitmapGenerator(100, ITExpRandomNumbers, 0);
% frequencyMonobitTest(ITExpRandomNumbers);
% runsTest(ITExpRandomNumbers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inverse Transform Sampling - Gauss
% parameters: number of numbers to generate, number of delay time, interval[0, 1]
% returns generated random numbers

% ITGaussRandomNumbers = inverseTransformGauss(1e4, 0, 0);
% bitmapGenerator(100, ITGaussRandomNumbers, 0);
% frequencyMonobitTest(ITGaussRandomNumbers);
% runsTest(ITGaussRandomNumbers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inverse Transform Sampling - Rayleigh
% parameters: number of numbers to generate, number of delay time, interval[0, 1]
% returns generated random numbers

% ITRayleighRandomNumbers = inverseTransformRayleigh(1e2, 0, 1);
% bitmapGenerator(100, ITRayleighRandomNumbers, 0);
% frequencyMonobitTest(ITRayleighRandomNumbers);
% runsTest(ITRayleighRandomNumbers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inverse Transform Sampling - Piecewise
% parameters: number of numbers to generate, number of delay time, interval[0, 1]
% returns generated random numbers

% ITPiecewiseRandomNumbers = inverseTransformPiecewise(1e2, 0, 1);
% bitmapGenerator(100, ITPiecewiseRandomNumbers, 0);
% frequencyMonobitTest(ITPiecewiseRandomNumbers);
% runsTest(ITPiecewiseRandomNumbers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bitmapGenerator(a, r, colorful)
    figure;
    matrice = zeros(a, a);
    n = 1;
    i = 1;
    j = 1;
    if ~colorful
       colormap([1 1 1; 0 0 0]);
       r = round(r);
    end
    
    while 1
        matrice(i, j) = r(n);
        n = n + 1;
        if n > length(r)
            n = 1;
        end
        i = i + 1;
        if i == a + 1
            i = 1;
            j = j + 1;
            if j == a + 1
                break;
            end
        end
    end
    imagesc(matrice);
    colorbar;
    axis square;
    title('Bitmap of Random Numbers');
end

function randomNumbers = linearCongruential(m, a, c, seed, number, isHist, isInterval)
    if ~exist('number','var')
        number = inf;
    end
    
    if seed >= m
        seed = mod(seed, m);
    end

    randomNumbers = [];
    n = 0;
    while ~sum(ismember(randomNumbers, seed))
        randomNumbers(end + 1) = seed;
        seed = mod(seed * a + c, m);
        n = n + 1;
        if n == number
            break;
            disp(['Period: ', num2str(length(randomNumbers))]);
        end
    end
    if isInterval
        randomNumbers = randomNumbers / max(randomNumbers); % 10^ceil(log10(max(randomNumbers)));
    end
    if isHist
        hist(randomNumbers);
    end
    disp(['Number of generated numbers: ', num2str(length(randomNumbers))]);
end
 
function X = inverseTransformExponential(n, time, isInterval)
    if ~exist('n','var')
        n = 10;
    end
    
    if ~exist('time','var')
        time = 0.01;
    end
    
    mu = 1.5;
    lambda = 1 / mu;
    x = 0:0.01:n;
    ypdf = lambda*exp(-lambda*x);       % f(x) % exppdf(x, mu);
    ycdf = 1 - exp(-lambda*x);          % F(x) % expcdf(x, mu);
    R = rand(1, n);                 % rand(1, 11);
    X = (-1/lambda)*log(1 - R);         % X = F^-1(R)
    
    figure;
    %hold on;
    plot(x, ypdf);
    axis([0 log(0.001/lambda)/-lambda 0 max(ypdf)]);
    xlabel('x');
    ylabel('f(x)');
    title('Probability Density Function (PDF)');
    
    figure
    plot(x, ycdf);
    xlabel('x');
    ylabel('F(x)');
    title('Cumulative Distrubution Function (CDF)');
    h = animatedline('Color','red','Marker','o','MarkerSize', 2,'MarkerFaceColor','r','LineStyle','none');

    i = 0;                   % X
    j = R(n);        % R
    k = n;
    
    if time ~= 0
        while 1
            if i >= X(k)
                addpoints(h, i, j);
                j = j - 1/200;
                if j <= 0
                    j = 0;
                    addpoints(h, i, j);
                    k = k - 1;
                    if k == 0
                        break;
                    end
                    j = R(k);
                    i = 0;
                end
            else
                addpoints(h, i, j);
                i = i + max(x)/200;
                if i >= X(k)
                    i = X(k);
                    addpoints(h, i, j);
                end
            end
            pause(time);
        end 
    end
    
    clf;
    hold on;
    plot(x, ycdf, 'Color', 'blue', 'LineWidth', 3);
    c = n;
    while c ~= 0
        line([0 X(c)], [R(c) R(c)], 'Color', 'red');
        line([X(c) X(c)], [R(c) 0], 'Color', 'red');
        c = c - 1;
    end
    plot(x, ycdf, 'Color', 'blue', 'LineWidth', 3);
    xlabel('x');
    ylabel('F(x)');
    title('Cumulative Distrubution Function (CDF)');
    axis([0 10 0 1]);
    
    figure;
    hist(R);
    title('Histogram of Uniform Distrubuted Random Numbers in the Interval (0, 1)');
    
    figure;
    histogram(X, 'Normalization', 'pdf');
    [f, xi] = ksdensity(X);
    hold on;
    plot(xi, f);
    title('Histogram of Exponential Distributed Generated Random Numbers');
    
    %X = X / 10;
    if isInterval
        X = X / 10^ceil(log10(max(X)));
    end
end

function X = inverseTransformGauss(n, time, isInterval)
    if ~exist('n','var')
        n = 10;
    end
    
    if ~exist('time','var')
        time = 0.01;
    end
    
    mean = 0;
    variance = 8;
    x = -10:0.01:10;
    ypdf = (1/sqrt(2*pi*variance))*exp((-(x-mean).^2)/(2*variance));
    ycdf = 0.5*(1+erf((x-mean)/(sqrt(variance*2))));                
    R = rand(1, n);                                                 
    X = mean + sqrt(variance*2)*(erfinv(2*R-1));                    
    
    figure;
    plot(x, ypdf);
    xlabel('x');
    ylabel('f(x)');
    title('Probability Density Function (PDF)');
    
    figure;
    plot(x, ycdf, 'Color', 'blue', 'LineWidth', 3);
    xlabel('x');
    ylabel('F(x)');
    title('Cumulative Distrubution Function (CDF)');
    h = animatedline('Color','red','Marker','o','MarkerSize', 2,'MarkerFaceColor','r','LineStyle','none');
    
    i = -10;                   % X
    j = R(n);        % R
    k = n;
  
    if time ~= 0
        while 1
            if i >= X(k)
                addpoints(h, i, j);
                j = j - 1/200;
                if j <= 0
                    j = 0;
                    addpoints(h, i, j);
                    k = k - 1;
                    if k == 0
                        break;
                    end
                    j = R(k);
                    i = -10;
                end
            else
                addpoints(h, i, j);
                i = i + max(x)/200;
                if i >= X(k)
                    i = X(k);
                    addpoints(h, i, j);
                end
            end
            pause(time);
        end 
    end
    clf;
    hold on;
    plot(x, ycdf, 'Color', 'blue', 'LineWidth', 3);
    c = n;
    while c ~= 0
        line([-10 X(c)], [R(c) R(c)], 'Color', 'red');
        line([X(c) X(c)], [R(c) 0], 'Color', 'red');
        c = c - 1;
    end
    plot(x, ycdf, 'Color', 'blue', 'LineWidth', 3);
    xlabel('x');
    ylabel('F(x)');
    title('Cumulative Distrubution Function (CDF)');
    axis([-10 10 0 1]);
    
    figure;
    hist(R);
    title('Histogram of Uniform Distrubuted Random Numbers in the Interval (0, 1)');
    
    figure;
    histogram(X, 'Normalization', 'pdf');
    [f, xi] = ksdensity(X);
    hold on;
    plot(xi, f);
    title('Histogram of Normal Distributed Generated Random Numbers');
    
    disp(['Number of generated numbers: ', num2str(n)]);
    disp(['Number of generated numbers greater than 0: ', num2str(sum(X > 0))]);
    disp(['Number of generated numbers less than 0: ', num2str(sum(X < 0))]);
    disp(['Number of generated numbers equal to 0: ', num2str(sum(X == 0))]);
    
    if isInterval
        X = X / 10^ceil(log10(max(X)));
    end
end

function X = inverseTransformRayleigh(n, time, isInterval)
    if ~exist('n','var')
        n = 10;
    end
    
    if ~exist('time','var')
        time = 0.01;
    end
    mode = 2;
    variance = mode^2*(4-pi)/2;
    x = 0:0.01:10;
    ypdf = exp(-x.^2/(2*mode.^2)).*x/mode.^2;       % f(x) % exppdf(x, mu);
    ycdf = 1-exp(-x.^2/(2*mode.^2));                       % F(x) % expcdf(x, mu);
    R = rand(1, n);                                                        % rand(1, n);
    X = mode*sqrt(-2*log(1-R));                           % X = F^-1(R)
    
    figure;
    plot(x, ypdf);
    xlabel('x');
    ylabel('f(x)');
    title('Probability Density Function (PDF)');
    
    figure;
    plot(x, ycdf, 'Color', 'blue', 'LineWidth', 3);
    xlabel('x');
    ylabel('F(x)');
    title('Cumulative Distrubution Function (CDF)');
    h = animatedline('Color','red','Marker','o','MarkerSize', 2,'MarkerFaceColor','r','LineStyle','none');
    
    i = 0;                   % X
    j = R(n);        % R
    k = n;
  
    if time ~= 0
        while 1
            if i >= X(k)
                addpoints(h, i, j);
                j = j - 1/200;
                if j <= 0
                    j = 0;
                    addpoints(h, i, j);
                    k = k - 1;
                    if k == 0
                        break;
                    end
                    j = R(k);
                    i = 0;
                end
            else
                addpoints(h, i, j);
                i = i + max(x)/200;
                if i >= X(k)
                    i = X(k);
                    addpoints(h, i, j);
                end
            end
            pause(time);
        end 
    end
    clf;
    hold on;
    plot(x, ycdf, 'Color', 'blue', 'LineWidth', 3);
    c = n;
    while c ~= 0
        line([-10 X(c)], [R(c) R(c)], 'Color', 'red', 'LineWidth', 0.1);
        line([X(c) X(c)], [R(c) 0], 'Color', 'red', 'LineWidth', 0.1);
        c = c - 1;
    end
    plot(x, ycdf, 'Color', 'blue', 'LineWidth', 3);
    xlabel('x');
    ylabel('F(x)');
    title('Cumulative Distrubution Function (CDF)');
    axis([0 10 0 1]);
    
    figure;
    hist(R);
    title('Histogram of Uniform Distrubuted Random Numbers in the Interval (0, 1)');
    
    figure;
    histogram(X, 'Normalization', 'pdf');
    [f, xi] = ksdensity(X);
    hold on;
    plot(xi, f);
    title('Histogram of Rayleigh Distributed Generated Random Numbers');
    
    disp(['Number of generated numbers: ', num2str(n)]);
    if isInterval
        X = X / 10^ceil(log10(max(X)));
    end
end

function X = inverseTransformPiecewise(n, time, isInterval)
    if ~exist('n','var')
        n = 10;
    end
    
    if ~exist('time','var')
        time = 0.01;
    end
    
    x = 0:0.01:1.5;
    ypdf1 = 2*x(1:100)/3;
    ypdf2 = (8*x(101:end)-6)/3;
    ypdf = [ypdf1, ypdf2];
    
    ycdf1 = x(1:100).^2/3;
    ycdf2 = 1/3 + (2*(2*x(101:end) - 1).*(x(101:end) - 1))/3;
    ycdf = [ycdf1, ycdf2];
    R = sort(rand(1, n));                                                        % rand(1, n);
    X1 = 3^(1/2).*R(R <= 1/3).^(1/2);
    X2 = (3*((16*R(R >= 1/3))./3 - 4/3).^(1/2))/8 + 3/4;
    X = [X1, X2];
    ix = randperm(n);
    X = X(ix);
    R = R(ix);
    
    figure;
    hold on;
    plot(x, ypdf);
    xlabel('x');
    ylabel('f(x)');
    title('Probability Density Function (PDF)');
    
    figure;
    hold on;
    plot(x, ycdf);
    xlabel('x');
    ylabel('F(x)');
    title('Cumulative Distrubution Function (CDF)');
    h = animatedline('Color','red','Marker','o','MarkerSize', 2,'MarkerFaceColor','r','LineStyle','none');
    
    i = 0;                   % X
    j = R(n);        % R
    k = n;
  
    if time ~= 0
        while 1
            if i >= X(k)
                addpoints(h, i, j);
                j = j - 1/200;
                if j <= 0
                    j = 0;
                    addpoints(h, i, j);
                    k = k - 1;
                    if k == 0
                        break;
                    end
                    j = R(k);
                    i = 0;
                end
            else
                addpoints(h, i, j);
                i = i + max(x)/200;
                if i >= X(k)
                    i = X(k);
                    addpoints(h, i, j);
                end
            end
            pause(time);
        end 
    end
    clf;
    hold on;
    plot(x, ycdf, 'Color', 'blue', 'LineWidth', 3);
    c = n;
    while c ~= 0
        line([-10 X(c)], [R(c) R(c)], 'Color', 'red', 'LineWidth', 0.1);
        line([X(c) X(c)], [R(c) 0], 'Color', 'red', 'LineWidth', 0.1);
        c = c - 1;
    end
    plot(x, ycdf, 'Color', 'blue', 'LineWidth', 3);
    xlabel('x');
    ylabel('F(x)');
    title('Cumulative Distrubution Function (CDF)');
    axis([0 1.5 0 1]);
    
    figure;
    hist(R);
    title('Histogram of uniform distrubuted random numbers in the interval (0, 1)');
    
    figure;
    histogram(X, 'Normalization', 'pdf');
    [f, xi] = ksdensity(X);
    hold on;
    plot(xi, f);
    title('Histogram of normal distributed generated random numbers');
    
    disp(['Number of generated numbers: ', num2str(n)]);
    
    if isInterval
        X = X / 10^ceil(log10(max(X)));
    end
end

function randomNumbers = middleSquare(seed, isHist, isInterval)
    if ~exist('seed','var')
        seed = input('Enter a random seed number: ');
    end

    disp(['Random seed: ', num2str(seed)]);
    seed = seed * 10^(~mod(ceil(log10(seed)), 2) == 0);
    n = ceil(log10(seed));

    randomNumbers = [];

    while ~sum(ismember(randomNumbers, seed))
        randomNumbers(end + 1) = seed;
        seed = seed^2;
        num = ceil(log10(seed));

        if num <= n+n/2
            if seed < 10^(n/2)
                seed = 0;
            else
                seed = floor(seed / 10^(n/2));
                
            end
        end

        for i = 1:n
            if num == n+i+n/2
                seed = floor( floor((seed - floor(seed / 10^(num-i)) * 10^(num-i))) / 10^(n/2));
                break;
            end
        end

    end
    period = length(randomNumbers);
    disp(['Period: ', num2str(period)]);
    if isInterval
        randomNumbers = randomNumbers / max(randomNumbers);%10^ceil(log10(max(randomNumbers)));
    end
    if isHist
        subplot(2, 1, 1);
        plot(randomNumbers);
        subplot(2, 1, 2);
        hist(randomNumbers);
    end
end

function frequencyMonobitTest(randomNumbers)
    disp('--------------------');
    disp('Frequency Monobit Test');
    format long;
    n = length(randomNumbers);
    numbers = round(randomNumbers);
    
    % a = num2str(numbers);
    % a(isspace(a)) = [];
    
    absSn = abs(n - 2 * sum(numbers));
    Sobs = absSn / sqrt(n);
    pValue = erfc(Sobs / sqrt(2));
    disp(['P-value: ', num2str(pValue)]);
    if pValue >= .01
        disp('P-value >= 0.01, accept the sequence as random.');
    else
        disp('P-value < 0.01, conclude that the sequence is non-random.');
    end
end

function runsTest(randomNumbers)
    disp('--------------------');
    disp('Runs Test');
    format long;
    n = length(randomNumbers);
    numbers = round(randomNumbers);
    pi = sum(numbers) / n;
    tau = 2 / sqrt(n);
    Vn = 0;
    if abs(pi - 0.5) >= tau
        disp('Because of pi - 0.5 < tau, then the Runs test are not be performed');
        pValue = 0;
    else
        for k = 1:n-1
            rk = numbers(k) ~= numbers(k + 1);
            Vn = Vn + rk;
        end
        Vn = Vn + 1;
        pValue = erfc(abs(Vn-2*n*pi*(1-pi))/(2*sqrt(2*n)*pi*(1-pi)));
    end
    
    disp(['P-value: ', num2str(pValue)]);
    if pValue >= 0.01
        disp('P-value >= 0.01, accept the sequence as random.');
    else
        disp('P-value < 0.01, conclude that the sequence is non-random.');
    end
end

