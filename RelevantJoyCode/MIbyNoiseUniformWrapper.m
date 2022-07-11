MIinfobyNoiseArray.noisevec = (0:.01:6);
for x = 1:7
    filename =[ 'Moth', num2str(x), '_MIdata.mat'];
    load(['D:\tenie_data\Desktop\Summer 2019\Moth Motor Timing\TobiasShare\SubmittedDataallmusclesAllareTzWsd/', filename])
    mothname = ['Moth', num2str(x)];
    fields = fieldnames(time_data);
    for fn=fields'
        fn;
        %# since +fn is a 1-by-1 cell array, you still need to index into it, unfortunately
        musclename = fn{1};
        musclename = musclename(1:end-7);
        MIvec = MutualInformationbynoiseuniform(filename, musclename);
        MIinfobyNoiseArray.(filename(1:5)).(musclename).MIvec =  MIvec;
        save('D:\tenie_data\Desktop\Summer 2019\Moth Motor Timing\TobiasShare\actual4pts\new\MIbyUniformnoise.mat', '-struct', 'MIinfobyNoiseArray')

        
    end
    
    
end