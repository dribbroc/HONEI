
#include <QtGui>

#include "glwidget.hh"
#include "window.hh"

Window::Window()
{
    glWidget = new GLWidget;
    menuwidget = new QLabel("bla");
/*
    m_xSlider = createSlider( 360 );
    m_ySlider = createSlider( 360 );
    m_flakeRecSlider = createSlider( 6 );
	m_backgroundColorSlider = createSlider( 100 );

    connect( m_xSlider, SIGNAL(valueChanged(int)), glWidget, SLOT(setXRotation(int)) );
    connect( glWidget, SIGNAL(xRotationChanged(int)), m_xSlider, SLOT(setValue(int)) );
    connect( m_ySlider, SIGNAL(valueChanged(int)), glWidget, SLOT(setYRotation(int)) );
    connect( glWidget, SIGNAL(yRotationChanged(int)), m_ySlider, SLOT(setValue(int)) );
    connect( m_flakeRecSlider, SIGNAL(valueChanged(int)), glWidget, SLOT(setFlakeRec(int)) );
    connect( m_backgroundColorSlider, SIGNAL(valueChanged(int)), glWidget, SLOT( setBackgroundColor(int)) );
*/
    QHBoxLayout *mainLayout = new QHBoxLayout;
    mainLayout->addWidget(menuwidget);
    mainLayout->addWidget(glWidget);
  /*  mainLayout->addWidget(m_xSlider);
    mainLayout->addWidget(m_ySlider);
    mainLayout->addWidget(m_flakeRecSlider);
    mainLayout->addWidget(m_backgroundColorSlider);*/
    setLayout(mainLayout);
/*
    m_xSlider->setValue( 15 );
    m_ySlider->setValue( 345 );
    m_flakeRecSlider->setValue( 2 );
    m_backgroundColorSlider->setValue( 10 );*/
    setWindowTitle( "Hello GL" );
}

