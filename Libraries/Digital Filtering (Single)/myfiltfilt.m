function y = myfiltfilt(filtObject,x)

y = flipud(filter(filtObject,flipud(filter(filtObject,x))));