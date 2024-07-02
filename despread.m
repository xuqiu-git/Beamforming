function X_d = despread(X, code)
    [M, N] = size(X);
    K = N/length(code);     % 信号被扩展的块数
    X_d = zeros(M,K);
    for i=1:M
        x = X(i,:);
        x = reshape(x, [length(code), K]);
        x_d = code*x;       % 按行依次解扩
        X_d(i,:) = x_d;
    end
end