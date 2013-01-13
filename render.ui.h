/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you want to add, delete, or rename functions or slots, use
** Qt Designer to update this file, preserving your code.
**
** You should not define a constructor or destructor in this file.
** Instead, write your code in functions called init() and destroy().
** These will automatically be called by the form's constructor and
** destructor.
*****************************************************************************/


void Render::SetOpacity()
{
    //std::cout << "newSlot" << std::endl;
	float s = opacity->value()/100.0;
	char buffer[10];
	sprintf_s(buffer,sizeof(buffer),"%.2f\n",s);
	printf("%s\n",buffer);
	opacityInfo->setText( tr( buffer ) );
}





