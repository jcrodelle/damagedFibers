% take in a spike train and injure it and output a new spike train:
function [output_train] = injureSpikeTrain(input_train, injury_type, injury_param)

if strcmp(injury_type,'normal_conduction') == 1
    output_train = input_train;
end

if strcmp(injury_type,'total_blocking') == 1
    output_train = 0*input_train;
end

if strcmp(injury_type,'delay_train') == 1
    delay = injury_param.delay;
    output_train = 0*input_train;
    output_train(delay:end) = input_train(1:end-delay+1);
end

if strcmp(injury_type,'anticipate_train') == 1
    anticipate = injury_param.anticipate;
    output_train = 0*input_train;
    output_train(1:end-anticipate+1) = input_train(anticipate:end);
end

if strcmp(injury_type,'increase_refract') == 1
    tau = injury_param.tau;
    
    [output_train] = increase_refract(input_train,tau);
end

if strcmp(injury_type,'interm_blocking') == 1
    freq = injury_param.interm_freq;
    output_train = 0*input_train;
    
    for n = 1:size(input_train,2)
        output_train(n) = 0.5*(sign(sin(freq*n))+1).*input_train(n);
    end
end

if strcmp(injury_type,'evoke_potentials') == 1
    output_train = input_train;
    evoked = injury_param.evoked;
    space = injury_param.evokedSpace;  
    evokeProb = injury_param.evokedProb; 
    
    for j = 1:length(input_train)
        if input_train(j) == 1 
            if rand(1) < evokeProb
                for k =1:evoked
                    if j+k*space<length(input_train)
                        output_train(j+k*space) = 1;
                    end
                end
            end
        end
    end
    
end
