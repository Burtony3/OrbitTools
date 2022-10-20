function [mu, sig, tvec] = btime(func, input, randmag, n)
    if nargin < 3; randmag = 1e-2; end
    if nargin < 4; n = 1e5; end
    tvec = zeros(n, 1);
    sz = size(input);

    textprogressbar('Profiling: ');
    for i = 1:n
        tic;
        func(input + randmag.*randn(sz));
        t = toc;
        tvec(i) = t;
        textprogressbar(i/n*100);
    end
    textprogressbar(' Done')
    mu = mean(tvec);
    sig = std(tvec);
end

