/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of HONEI. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#include <QApplication>

#include <clients/qt/window.hh>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    Window window;
    window.show();
    return app.exec();
}

