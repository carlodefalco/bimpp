

clear all
close all

load matdata
maxts = 100;
for jj=0:maxts

	fprintf(1,'frame %d\n',jj);

	filename = num2mstr(sprintf('step_%5.5d',jj));
	vecname = num2mstr(sprintf('step%5.5d',jj));
	eval(filename);
	
	pdeplot(p,e,t,'xydata',eval(vecname),'zdata',eval(vecname),'colormap','jet');
	%set(gcf,'visible','off')
	caxis([ 0 1])
	axis([0 1 0 1 0 1])
	colorbar

	box on
	grid on

	print('-dpng',['frame' num2str(jj)])
	pause(1)
        close all
end