filename = 'testAnimated.gif';
colormap jet

for n = 1:100

    h = figure('position', [200 200 600 600]);
    isosurf(X,Y,Z,abs(Adata(:,:,:,n)))
    %drawnow 
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
    close(h)
  end

