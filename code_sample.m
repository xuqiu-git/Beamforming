function sample_code = code_sample(original_code,fc,fs,Ns)
% Ns��������,fs����Ƶ��,fc�����Ƶ��
    index = (1 : Ns)/(fs/fc);
    n = floor(mod(index,fc*0.001)) + 1;
    sample_code=original_code(n);
end