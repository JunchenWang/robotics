function blkStruct = slblocks
    Browser.Library = 'my_library';
    Browser.Name = 'MyLibrary';
    Browser.IsFlat = 0;

    blkStruct.Browser = Browser;
end