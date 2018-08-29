% PJC wrote this
% plot3(collection_of_chains(:,1,i),collection_of_chains(:,2,i),collection_of_chains(:,3,i))

figure;
[r,c,p]=size(collection_of_chains);

% for every crystal
ps=collection_of_chains;
for i=1:p
    h=patch(ps(1:6,1,i),ps(1:6,2,i),ps(1:6,3,i),'b');
    set(h,'facealpha',0.5)
    h=patch(ps([1:6]+6,1,i),ps([1:6]+6,2,i),ps([1:6]+6,3,i),'b');hold on;
    set(h,'facealpha',0.5)
    for j=1:5
        h=patch([ps([1:2]+(j-1),1,i); ps([2:-1:1]+6+(j-1),1,i)],...
            [ps([1:2]+(j-1),2,i); ps([2:-1:1]+6+(j-1),2,i)],...
            [ps([1:2]+(j-1),3,i); ps([2:-1:1]+6+(j-1),3,i)],'b');
        set(h,'facealpha',0.5)
    end
    h=patch([ps([6 1],1,i); ps([7 12],1,i)],...
        [ps([6 1],2,i); ps([7 12],2,i)],...
        [ps([6 1],3,i); ps([7 12],3,i)],'b');
    set(h,'facealpha',0.5)
end
set(gcf,'renderer','openGL')
lighting phong
light
axis equal
% set(gca,'visible','off')