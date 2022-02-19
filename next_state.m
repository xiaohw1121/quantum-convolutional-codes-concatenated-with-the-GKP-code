function [state_next,memory_contents]=next_state(state,input,L)
state = dec2bin(state,L-1);
state_next = zeros(1,length(state));
state_next(2:length(state)) = state(1:length(state)-1);
if input == 0
    state_next(1) = '0';
    else
    state_next(1) = '1';
end
memory_contents = [state_next,state(length(state))];
state_next = state_next - 48;
memory_contents = memory_contents -48;
y = 0;
for i = 1:length(state_next)
    y = y+state_next(length(state_next)-i+1).*2.^(i-1);
end
state_next = y;

end

