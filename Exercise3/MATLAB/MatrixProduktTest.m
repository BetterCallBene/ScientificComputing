
h = 1.0/4.0;
N = round(1.0/h + 1.0);

iK = (N-2)^2;

A= zeros(iK, iK);

curElm = 1;
lastElm = (N-2);

for i = 1:1:lastElm
    North =(i ~= 1);
    South = (i ~= lastElm);
    for j= 1:1:lastElm
		
        West = (j ~= 1);
        East = (j ~= lastElm);

        curElm = (i-1) * lastElm + j;	
        A(curElm, curElm) = 4;
        if West == true
            A(curElm,  curElm-1) = -1;
        end
        if East == true
			A(curElm,  curElm+1) = -1;
        end
        if South == true
            A(curElm, curElm + lastElm) = -1;
        end
        if North == true
			A(curElm, curElm - lastElm) = -1;	
        end
    end
end
            