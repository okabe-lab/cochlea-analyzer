function cellLossMat = SimulatePositionalEffect(COCHLEAR_SIZE, LOSS_RATIO...
    , powerInd, gaussSigma, TRIAL_NO)
%*****************************************************************************
%This function uses "COCHLEAR_SIZE", numbers of rows and columns of virtual
%outer hair cells, "LOSS_RATIO", measured ratio of cell loss, "powerInd",
%parameter of power-law distribution, "gaussSigma", standard deviation of
%Gaussian smooting kernel, and "TRIAL_NO". The output is "cellLossMat",
%result of simulations indicating positions of cell loss.
%*****************************************************************************

% Simulate cell loss with "Positional Effect", which is created on specified 
% parameters of power-low distribuion and Gaussian filtering.
cellLossMat = zeros(COCHLEAR_SIZE(1), COCHLEAR_SIZE(2), TRIAL_NO);
LostCell = floor(COCHLEAR_SIZE(1)*COCHLEAR_SIZE(2)*LOSS_RATIO);

for i=1:1000
    Power= 0.1*powerInd*rand(COCHLEAR_SIZE(2),1).^-(1-0.1*powerInd);
    GPower =imgaussfilt(Power,gaussSigma);
    GPower = GPower/sum(GPower);   
    SumGPower = cumsum(GPower);
    
    judge = 0;
    Ncount = 0;
    while judge == 0
        positionY = (ceil(COCHLEAR_SIZE(1)*rand));
        positionX = find(SumGPower > rand, 1);
        if cellLossMat(positionY, positionX, i) ==0
            cellLossMat(positionY, positionX, i) =1;
            Ncount = Ncount +1;
        end
        if Ncount == LostCell
            judge= 1;
        end
    end
end