
function error=segmentationError(groundTrueLabel,label)
groundTrueLabel=groundTrueLabel+ones(1,length(groundTrueLabel));
label=label+ones(1,length(label));
number=length(unique(groundTrueLabel));
flag=zeros(1,max(groundTrueLabel));
% label=groundTrueLabel;
n=length(unique(label));
rightnumber=0;
for i=1:max(label)
    A(i)=sum(label == i);
end

[ss,sortnumber]=sort(A, 'descend');

for i=1:max(label)
    each_cluster=(label == sortnumber(i));
    maxx = 0;
    kk=-1;
    for j=min(groundTrueLabel):max(groundTrueLabel)
        if flag(j)==1 
            continue
        end
        True_cluster=(groundTrueLabel == j);
        k=length(find(each_cluster+True_cluster>1));
        if k> maxx
            maxx = k; 
            kk=j;
        end
    end
    
    if maxx>0
        flag(kk)=1;
        rightnumber=rightnumber+maxx;
    end  
end

error=(length(label)-rightnumber)/length(label);

