function [t1,t2] = Switch_selector(signal_num)
     if signal_num<=1
        t1=1;
        t2=1;
    elseif signal_num<=2
        t1=0;
        t2=1;
     else
        t1=0;
        t2=0;
    end

end

