Summary: Fortran 90/95 bindings for shapelib
Name: shapelib-fortran
Version: 1.3
Release: 1
License: LGPL
Group: Deelopment/Libraries
URL: http://shapelib.maptools.org/
Source0: %{name}-%{version}.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
BuildRequires: shapelib-devel

%package devel

Group: Deelopment/Libraries
Summary: Static libraries and modules for shapelib-fortran.

%description

Bindings for calling shapelib API from fortran 90/95.

%description devel

Bindings for calling shapelib API from fortran 90/95 devel package.

%prep
%setup -q

%build
%configure
make

%install
rm -rf %{buildroot}
make DESTDIR=%{buildroot} install

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root,-)
%{_libdir}/*.la
%{_libdir}/*.so*
##%{_docdir}/%{name}-%{version}/*
%doc LICENSE shapelib-fortran.html

%files devel
%{_includedir}/*.mod
%{_libdir}/*.a

%changelog
* Wed Jun  1 2011 Davide Cesari <cesari@malina.metarpa> - 1.3-1
- New version

* Wed Oct  1 2008 Davide Cesari <cesari@malina.metarpa> - 1.2-1
- New version

* Thu Apr 17 2008 Davide Cesari <cesari@malina.metarpa> - 1.1-1
- New version

* Mon Dec 10 2007 Davide Cesari <cesari@malina.metarpa> - 1.0-3
- Fixed makeinstall for Fedora 8.

* Mon Nov 19 2007 Davide Cesari <cesari@malina.metarpa> - 1.0-2
- Added devel package and dependencies.

* Wed Apr 11 2007 Davide Cesari <cesari@malina.metarpa> 1.0-1
- Initial build.

