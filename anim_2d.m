filename = 'testAnimated.gif';
h = figure('position', [200 200 600 600]);
colormap jet

for n = 1:101
    % Draw plot for y = x.^n
    
    surf(X,Y, abs(Adata(:,:,n)))
    view(0,90), shading interp, axis tight
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
  end

