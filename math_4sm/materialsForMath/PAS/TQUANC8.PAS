program Tquanc8;
{$N-}
uses FMM, crt;
Var
   a, b, relerr, abserr,
   result, errest, flag, k : float;
   nofun,i                 : longint;
{$F+}
Function fun(y : float) : float;
begin
  if y=0.0 then fun:=k
           else fun:=sin(k*y)/y;
end;
{$F-}
begin
 clrscr;
 writeln('  ??????p?p????? ?p??p???? ??? QUANC8');
 writeln;
 writeln(' ?㤥? ??⥣p?p????? sin(k * x) / x ');
 writeln(' ???????, ??????????, ?।??? ??⥣??஢????:');
 write('              ??????,  a:= ');readln(a);
 write('              ???孨?, b:= ');readln(b);
 relerr:=1.0e-6;
 abserr:=0.0;
 k:=1;
 writeln;
 writeln('  K        Result          Errest         flag       Nofun');
 writeln;
 while k < 101.1 do begin
            quanc8(@fun,a,b,abserr,relerr,result,errest,nofun,flag);
            writeln(k:3:0, '  ',result:14:10,'  ',errest,'   ',flag:5:2,
             '       ',nofun:4);
            k:=k+10;
                  end;
end.
