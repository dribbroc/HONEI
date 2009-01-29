
#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QLabel>
#include "glwidget.hh"

class QSlider;
class GLWidget;

class Window : public QWidget
{
    Q_OBJECT

public:
    Window();

private:
  //  QSlider *createSlider( unsigned int range );

    GLWidget *glWidget;
    QLabel *menuwidget;
/*    QSlider *m_xSlider;
    QSlider *m_ySlider;
    QSlider *m_flakeRecSlider;
    QSlider *m_backgroundColorSlider;*/
};

#endif
