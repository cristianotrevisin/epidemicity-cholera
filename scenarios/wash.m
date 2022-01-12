 load ../data/wash date_list_cati cati_day
 
 a = datenum(2013,06,01);
b = datenum(2017,07,31);

figure()
hold on
for i=1:10
    bar(date_list_cati(2:end),cati_day(:,i),'stacked')
end
set(gca,'Xlim',[a b])
    datetick('x','mmm-yy','keeplimits','keepticks')
    box off
    legend("Artibonite", "Centre", "Grande Anse", "Nippes", "Nord-Est", "Nord-Ouest", "Ouest", "Sud", "Sud-Est", "Location","bestoutside")
    legend boxoff
    box off
    f.Units='points';
    f.Position=[0 0 400 400];