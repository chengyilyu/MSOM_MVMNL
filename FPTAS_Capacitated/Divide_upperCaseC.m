function FinalDivisionMatix = Divide_upperCaseC(numberofGroup,upperCaseC,MaxLowerCasec,MinLowerCasec) 
    if (MinLowerCasec*numberofGroup > upperCaseC) || MaxLowerCasec*numberofGroup < upperCaseC
        DivisionMatix = [];
    else
        p = MaxLowerCasec - MinLowerCasec;
        n = upperCaseC - numberofGroup*MinLowerCasec;
        s = {[0,1]};
        for i = 2 : n     
            k = size(s,2);
            w = {};
            for j = 1 : k
                if (s{j}(end) < s{j}(end-1)) && (s{j}(end) < p)
                    w{end+1} = s{j};
                    w{end}(end) = w{end}(end)+1;
                end
                if size(s{j},2) < numberofGroup+1
                    w{end+1} = s{j};
                    w{end}(end+1) = 1;
                end
            end
            if(i <= p)
                w{end+1} = [0,i];
            end
            s=w;
        end
        if numberofGroup*MinLowerCasec == upperCaseC    %(1)
            s={[0]};
        end
        DivisionMatix = ones(size(s,2),numberofGroup)*MinLowerCasec;
        for i=1:size(s,2)
            DivisionMatix(i,1:size(s{i},2)-1)=s{i}(2:end) + MinLowerCasec;
        end
        
        FinalDivisionMatix = [];
        for index = 1 : size(DivisionMatix,1)
            tempMatrix = [];
            tempMatrix = perms(DivisionMatix(index,:));
            FinalDivisionMatix = [FinalDivisionMatix; tempMatrix];
        end
        FinalDivisionMatix = unique(FinalDivisionMatix,"rows","sorted");
    end
end