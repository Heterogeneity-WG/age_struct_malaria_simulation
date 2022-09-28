load('scatter_xyz.mat')

figure()
subplot(1,2,1)
imagesc(xxs([1,end],1),yys(1,[1,end]),zz)
xlabel('x'); ylabel('y')
set(gca,'YDir','normal')
title('imagesc with x and y - two element vector')
subplot(1,2,2)
imagesc(xxs(:,1),yys(1,:),zz)
xlabel('x'); ylabel('y')
set(gca,'YDir','normal')
title('imagesc with x and y full vector')

%%
figure()
subplot(2,2,1)
imagesc(xxs(:,1),yys(1,:),zz)
xlabel('x'); ylabel('y')
set(gca,'YDir','normal')
title('imagesc with x and y specified')

subplot(2,2,2)
surf(xxs,yys,zz)
shading interp
xlabel('x'); ylabel('y')
title('surf with x and y specified')
view(3)

subplot(2,2,3)
surf(zz)
title('flattened surf without x and y')
shading interp
view(0,90)

subplot(2,2,4)
surf(zz)
title('surf without x and y')
shading interp
view(3)