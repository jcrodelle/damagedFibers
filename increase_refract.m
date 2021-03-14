function [output_train] = increase_refract(input_train,tau)

output_train = 0*input_train;
output_train(1) = input_train(1);

for j =1:tau
    for n = j+1: length(input_train)
        my_arg = input_train(n) + output_train(n-j) -2;
        
        if my_arg==0
            my_kron = 1;
        else
            my_kron = 0;
        end
        
        output_train(n) = input_train(n) - my_kron;
    end
    input_train = output_train;
end

end